using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers.MenkBordas;
using ISAAR.MSolve.XFEM.Utilities;

//TODO: This needs to be split between standard and enriched. Otherwise all enriched stuff may be uninitialized or in a 
//      previous state.
namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    class XClusterDofOrderer: IDofOrderer
    {
        private readonly XCluster2D cluster;
        private readonly DofTable<DisplacementDof> constrainedDofs;
        private readonly DofTable<DisplacementDof> standardDofs;

        private XClusterDofOrderer(XCluster2D cluster, int numConstrainedDofs,
            DofTable<DisplacementDof> constrainedDofs, int numStandardDofs, DofTable<DisplacementDof> standardDofs)
        {
            this.cluster = cluster;
            this.NumConstrainedDofs = numConstrainedDofs;
            this.constrainedDofs = constrainedDofs;
            this.NumStandardDofs = numStandardDofs;
            this.standardDofs = standardDofs;
        }

        public int NumConstrainedDofs { get; }
        public int NumEnrichedDofs
        {
            get
            {
                int count = 0;
                foreach (var subdomain in cluster.Subdomains)
                {
                    if (subdomain.DofOrderer != null) count += subdomain.DofOrderer.NumEnrichedDofs;
                }
                return count;
            }
        }
        public int NumStandardDofs { get; }


        /// <summary>
        /// Only standard free and standard constrained dofs.
        /// </summary>
        /// <param name="model"></param>
        /// <param name="cluster"></param>
        /// <returns></returns>
        public static XClusterDofOrderer CreateNodeMajor(Model2D model, XCluster2D cluster)
        {
            (int numConstrainedDofs, DofTable<DisplacementDof> constrainedDofs) = OrderConstrainedDofs(model.Constraints);
            (int numStandardDofs, DofTable<DisplacementDof> standardDofs) = OrderStandardDofs(model);

            return new XClusterDofOrderer(cluster, numConstrainedDofs, constrainedDofs, numStandardDofs, standardDofs);
        }

        public IDofOrderer DeepCopy()
        {
            //TODO: for this to work, there must not be stored references to the cluster or dubdomains, but their DofOrderers
            //      directly.
            throw new NotImplementedException();
        }

        /// <summary>
        /// TODO: Modify this method to extract any kind of vector, not only displacement, which means different 
        /// handling of constrained dofs, if constrained dofs are defined in the first place.
        /// </summary>
        /// <param name="element"></param>
        /// <param name="globalFreeVector"></param>
        /// <param name="globalConstrainedVector"></param>
        /// <returns></returns>
        public Vector ExtractDisplacementVectorOfElementFromGlobal(XContinuumElement2D element,
            Vector globalFreeVector, Vector globalConstrainedVector)
        {
            DofTable<DisplacementDof> elementDofs = element.GetStandardDofs();
            double[] elementVector = new double[elementDofs.EntryCount];
            foreach (Tuple<XNode2D, DisplacementDof, int> entry in elementDofs)
            {
                bool isStandard = this.standardDofs.TryGetValue(entry.Item1, entry.Item2, out int globalStandardDof);
                if (isStandard) elementVector[entry.Item3] = globalFreeVector[globalStandardDof];
                else
                {
                    int globalConstrainedDof = this.constrainedDofs[entry.Item1, entry.Item2];
                    elementVector[entry.Item3] = globalConstrainedVector[globalConstrainedDof];
                }
            }
            return Vector.CreateFromArray(elementVector);
        }

        /// <summary>
        /// </summary>
        /// <param name="element"></param>
        /// <param name="globalFreeVector">Both the free standard and enriched dofs.</param>
        /// <returns></returns>
        public Vector ExtractEnrichedDisplacementsOfElementFromGlobal(XContinuumElement2D element, Vector globalFreeVector)
        {
            XSubdomain2D subdomain = cluster.FindSubdomainOfElement(element);
            return subdomain.DofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(cluster, element, globalFreeVector);
        }

        public ITable<XNode2D, EnrichedDof, double> GatherEnrichedNodalDisplacements(Model2D model, Vector solution)
        {
            throw new InvalidOperationException("This method does not make sense for this dofOrderer. Refactor them.");
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="solution"></param>
        /// <returns>A nodesCount x 2 array, where each row stores the x and y displacements of that node</returns>
        public double[,] GatherNodalDisplacements(Model2D model, Vector solution)
        {
            double[,] result = new double[model.Nodes.Count, 2];
            for (int i = 0; i < model.Nodes.Count; ++i)
            {
                XNode2D node = model.Nodes[i];

                bool isXStandard = standardDofs.TryGetValue(node, DisplacementDof.X, out int globalStandardDofX);
                if (isXStandard) result[i, 0] = solution[globalStandardDofX];
                else result[i, 0] = model.Constraints[node, DisplacementDof.X];

                bool isYStandard = standardDofs.TryGetValue(node, DisplacementDof.Y, out int globalStandardDofY);
                if (isYStandard) result[i, 1] = solution[globalStandardDofY];
                else result[i, 1] = model.Constraints[node, DisplacementDof.Y];
            }
            return result;
        }

        public int GetConstrainedDofOf(XNode2D node, DisplacementDof dofType)
        {
            return constrainedDofs[node, dofType];
        }

        public IEnumerable<int> GetConstrainedDofsOf(XNode2D node)
        {
            // Perhaps it would be more efficient for the client to traverse all (dofType, dofIdx) pairs 
            // than assembling the dofIdx collection beforehand.
            try
            {
                return constrainedDofs.GetValuesOfRow(node);
            }
            catch (KeyNotFoundException)
            {
                return new int[] { }; // TODO: It would be better if this method is not called for a non enriched node.
            }
        }

        public int GetEnrichedDofOf(XNode2D node, EnrichedDof dofType)
        {
            throw new InvalidOperationException("This method does not make sense for this dofOrderer. Refactor them.");
        }

        public IEnumerable<int> GetEnrichedDofsOf(XNode2D node)
        {
            throw new InvalidOperationException("This method does not make sense for this dofOrderer. Refactor them.");
        }

        public List<int> GetEnrichedDofsOf(XContinuumElement2D element)
        {
            throw new InvalidOperationException("This method does not make sense for this dofOrderer. Refactor them.");
        }

        public List<int> GetConstrainedDofsOf(XContinuumElement2D element)
        {
            var globalDofs = new List<int>();
            foreach (var nodeDofLocal in element.GetStandardDofs())
            {
                bool isConstrained = constrainedDofs.TryGetValue(nodeDofLocal.Item1, nodeDofLocal.Item2, out int globalDof);
                if (isConstrained) globalDofs.Add(globalDof);
            }
            return globalDofs;
        }

        public int GetStandardDofOf(XNode2D node, DisplacementDof dofType)
        {
            return standardDofs[node, dofType];
        }

        public IEnumerable<int> GetStandardDofsOf(XNode2D node)
        {
            // Perhaps it would be more efficient for the client to traverse all (dofType, dofIdx) pairs 
            // than assembling the dofIdx collection beforehand.
            try
            {
                return standardDofs.GetValuesOfRow(node);
            }
            catch (KeyNotFoundException)
            {
                return new int[] { }; // TODO: It would be better if this method is not called for a non enriched node.
            }
        }

        public List<int> GetStandardDofsOf(XContinuumElement2D element)
        {
            var globalDofs = new List<int>(2 * element.Nodes.Count);
            foreach (var nodeDofLocal in element.GetStandardDofs())
            {
                bool isStandard = standardDofs.TryGetValue(nodeDofLocal.Item1, nodeDofLocal.Item2, out int globalDof);
                if (isStandard) globalDofs.Add(globalDof);
            }
            return globalDofs;
        }

        /// <summary>
        /// Index i = element local dof. Dictionary[i] = global dof.
        /// </summary>
        /// <param name="element"></param>
        /// <returns></returns>
        public IReadOnlyDictionary<int, int> MatchElementToGlobalEnrichedDofsOf(XContinuumElement2D element) //TODO: this can be an array. 
        {
            XSubdomain2D subdomain = cluster.FindSubdomainOfElement(element);
            return subdomain.DofOrderer.MatchElementToGlobalEnrichedDofs(element);
        }

        /// <summary>
        /// Index i = element local dof. Dictionary[i] = global dof. 
        /// </summary>
        /// <param name="element"></param>
        /// <returns></returns>
        public void MatchElementToGlobalStandardDofsOf(XContinuumElement2D element,
            out IReadOnlyDictionary<int, int> elementToGlobalStandardDofs,
            out IReadOnlyDictionary<int, int> elementToGlobalConstrainedDofs)
        {
            DofTable<DisplacementDof> elementDofs = element.GetStandardDofs();
            var globalStandardDofs = new Dictionary<int, int>();
            var globalConstrainedDofs = new Dictionary<int, int>();
            foreach (Tuple<XNode2D, DisplacementDof, int> entry in elementDofs)
            {
                bool isStandard = this.standardDofs.TryGetValue(entry.Item1, entry.Item2, out int standardGlobalDof);
                if (isStandard) globalStandardDofs[entry.Item3] = standardGlobalDof;
                else globalConstrainedDofs[entry.Item3] = this.constrainedDofs[entry.Item1, entry.Item2];
            }
            elementToGlobalStandardDofs = globalStandardDofs;
            elementToGlobalConstrainedDofs = globalConstrainedDofs;
        }

        public void OrderSubdomainDofs(ISet<XSubdomain2D> enrichedSubdomains, ISet<XSubdomain2D> modifiedSubdomains,
            ICrackDescription crack)
        {
            int numTotalDofs = NumStandardDofs;
            foreach (var subdomain in enrichedSubdomains)
            {
                if (modifiedSubdomains.Contains(subdomain))
                {
                    subdomain.DofOrderer = XSubdomainDofOrderer.CreateNodeMajor(crack, subdomain);
                }
                subdomain.DofOrderer.FirstGlobalDofIndex = numTotalDofs;
                numTotalDofs += subdomain.DofOrderer.NumEnrichedDofs;
            }
        }

        public void ReorderStandardDofs(IReadOnlyList<int> permutation, bool oldToNew)
        {
            standardDofs.Reorder(permutation, oldToNew);
        }

        /// <summary>
        /// Renumbers the dof indices according th the given permutation vector and direction. 
        /// If (<paramref name="oldToNew"/> == true), then newIndex[dof] = <paramref name="permutation"/>[oldIndex[dof]].
        /// Else oldIndex[dof] = <paramref name="permutation"/>[nwIndex[dof]]
        /// </summary>
        /// <param name="permutation">The permutation vector.</param>
        /// <param name="oldToNew">The direction it should be applied to.</param>
        public void ReorderUnconstrainedDofs(IReadOnlyList<int> permutation, bool oldToNew)
        {
            // Reorder standard free dofs
            // Ask subdomains to reorder their enriched dofs
            throw new NotImplementedException();
        }

        public void WriteToConsole()
        {
            Console.WriteLine("Standard free dofs for the whole domain:");
            Console.WriteLine(standardDofs);
            Console.WriteLine("Standard constrained dofs for the whole domain:");
            Console.WriteLine(constrainedDofs);
            for (int s = 0; s < cluster.Subdomains.Count; ++s)
            {
                Console.WriteLine($"Subdomain {s}:");
                cluster.Subdomains[s].DofOrderer.WriteToConsole();
            }
        }

        private static (int numConstrainedDofs, DofTable<DisplacementDof> constrainedDofs) OrderConstrainedDofs(
            ITable<XNode2D, DisplacementDof, double> constraints)
        {
            var constrainedDofs = new DofTable<DisplacementDof>();
            int counter = 0;
            foreach (Tuple<XNode2D, DisplacementDof, double> entry in constraints)
            {
                constrainedDofs[entry.Item1, entry.Item2] = counter++;
            }
            return (counter, constrainedDofs);
        }

        private static (int numStandardDofs, DofTable<DisplacementDof> standardDofs) OrderStandardDofs(Model2D model)
        {
            ITable<XNode2D, DisplacementDof, double> constraints = model.Constraints;
            var standardDofs = new DofTable<DisplacementDof>();
            int counter = 0;
            foreach (var node in model.Nodes)
            {
                if (!constraints.Contains(node, DisplacementDof.X)) standardDofs[node, DisplacementDof.X] = counter++;
                if (!constraints.Contains(node, DisplacementDof.Y)) standardDofs[node, DisplacementDof.Y] = counter++;
            }
            return (counter, standardDofs);
        }
    }
}
