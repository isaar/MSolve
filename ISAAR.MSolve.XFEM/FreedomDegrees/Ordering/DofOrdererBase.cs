using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Utilities;

//TODO: I should use interchangeble builders, rather than subclasses.
//TODO: constrained dofs should be ignored (set to -1) and the effect of constraints should be done via equivalent forces.
namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    class DofOrdererBase: IDofOrderer
    {
        protected readonly DofTable<DisplacementDof> constrainedDofs;
        protected readonly DofTable<EnrichedDof> enrichedDofs;
        protected readonly DofTable<DisplacementDof> freeDofs;

        public int ConstrainedDofsCount { get; }
        public int EnrichedDofsCount { get; }
        public int FreeDofsCount { get; }

        protected DofOrdererBase(int constrainedDofsCount, DofTable<DisplacementDof> constrainedDofs, 
            int enrichedDofsCount, DofTable<EnrichedDof> enrichedDofs,
            int freeDofsCount, DofTable<DisplacementDof> freeDofs)
        {
            this.ConstrainedDofsCount = constrainedDofsCount;
            this.constrainedDofs = constrainedDofs;
            this.EnrichedDofsCount = enrichedDofsCount;
            this.enrichedDofs = enrichedDofs;
            this.FreeDofsCount = freeDofsCount;
            this.freeDofs = freeDofs;
        }

        public IDofOrderer DeepCopy()
        {
            return new DofOrdererBase(ConstrainedDofsCount, constrainedDofs.DeepCopy(),
                EnrichedDofsCount, enrichedDofs.DeepCopy(), FreeDofsCount, freeDofs.DeepCopy());
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
            ITable<XNode2D, DisplacementDof, int> elementDofs = element.GetStandardDofs();
            double[] elementVector = new double[elementDofs.EntryCount];
            foreach (Tuple<XNode2D, DisplacementDof, int> entry in elementDofs)
            {
                bool isFree = this.freeDofs.TryGetValue(entry.Item1, entry.Item2, out int globalFreeDof);
                if (isFree) elementVector[entry.Item3] = globalFreeVector[globalFreeDof];
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
        /// <param name="globalFreeVector">Both the free and enriched dofs.</param>
        /// <returns></returns>
        public Vector ExtractEnrichedDisplacementsOfElementFromGlobal(XContinuumElement2D element,
            Vector globalFreeVector)
        {
            ITable<XNode2D, EnrichedDof, int> elementDofs = element.GetEnrichedDofs();
            double[] elementVector = new double[elementDofs.EntryCount];
            foreach (Tuple<XNode2D, EnrichedDof, int> entry in elementDofs)
            {
                int globalFreeDof = enrichedDofs[entry.Item1, entry.Item2];
                elementVector[entry.Item3] = globalFreeVector[globalFreeDof];
            }
            return Vector.CreateFromArray(elementVector);
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

        public int GetEnrichedDofOf(XNode2D node, EnrichedDof dofType)
        {
            return enrichedDofs[node, dofType];
        }

        public IEnumerable<int> GetEnrichedDofsOf(XNode2D node)
        {
            try
            {
                return enrichedDofs.GetValuesOfRow(node);
            }
            catch (KeyNotFoundException)
            {
                return new int[] { }; // TODO: It would be better if this method is not called for a non enriched node.
            }
        }

        public List<int> GetEnrichedDofsOf(XContinuumElement2D element)
        {
            var globalDofs = new List<int>();
            foreach (var nodeDofLocal in element.GetEnrichedDofs())
            {
                globalDofs.Add(enrichedDofs[nodeDofLocal.Item1, nodeDofLocal.Item2]); // It must be included.
            }
            return globalDofs;
        }

        // Would it be faster to return the Dictionary<StandardDofType, int> for consecutive accesses of the dofs of this node? 
        // Dictionary's Read operation is supposed to be O(1), but still...
        public int GetFreeDofOf(XNode2D node, DisplacementDof dofType)
        {
            return freeDofs[node, dofType];
        }

        public IEnumerable<int> GetFreeDofsOf(XNode2D node)
        {
            // Perhaps it would be more efficient for the client to traverse all (dofType, dofIdx) pairs 
            // than assembling the dofIdx collection beforehand.
            try
            {
                return freeDofs.GetValuesOfRow(node);
            }
            catch (KeyNotFoundException)
            {
                return new int[] { }; // TODO: It would be better if this method is not called for a non enriched node.
            }
        }

        public List<int> GetFreeDofsOf(XContinuumElement2D element)
        {
            var globalDofs = new List<int>(2 * element.Nodes.Count);
            foreach (var nodeDofLocal in element.GetStandardDofs())
            {
                bool isFree = freeDofs.TryGetValue(nodeDofLocal.Item1, nodeDofLocal.Item2, out int globalDof);
                if (isFree) globalDofs.Add(globalDof);
            }
            return globalDofs;
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

                bool isXFree = freeDofs.TryGetValue(node, DisplacementDof.X, out int globalFreeDofX);
                if (isXFree) result[i, 0] = solution[globalFreeDofX];
                else result[i, 0] = model.Constraints[node, DisplacementDof.X];

                bool isYFree = freeDofs.TryGetValue(node, DisplacementDof.Y, out int globalFreeDofY);
                if (isYFree) result[i, 1] = solution[globalFreeDofY];
                else result[i, 1] = model.Constraints[node, DisplacementDof.Y];
            }
            return result;
        }

        public ITable<XNode2D, EnrichedDof, double> GatherEnrichedNodalDisplacements(Model2D model,
            Vector solution)
        {
            var table = new Table<XNode2D, EnrichedDof, double>();
            foreach (var row in enrichedDofs)
            {
                table[row.Item1, row.Item2] = solution[row.Item3];
            }
            return table;
        }

        /// <summary>
        /// Index i = element local dof. Dictionary[i] = global dof. 
        /// </summary>
        /// <param name="element"></param>
        /// <returns></returns>
        public void MatchElementToGlobalStandardDofsOf(XContinuumElement2D element,
            out IReadOnlyDictionary<int, int> elementToGlobalFreeDofs,
            out IReadOnlyDictionary<int, int> elementToGlobalConstrainedDofs)
        {
            ITable<XNode2D, DisplacementDof, int> elementDofs = element.GetStandardDofs();
            var globalFreeDofs = new Dictionary<int, int>();
            var globalConstrainedDofs = new Dictionary<int, int>();
            foreach (Tuple<XNode2D, DisplacementDof, int> entry in elementDofs)
            {
                bool isFree = this.freeDofs.TryGetValue(entry.Item1, entry.Item2, out int freeGlobalDof);
                if (isFree) globalFreeDofs[entry.Item3] = freeGlobalDof;
                else globalConstrainedDofs[entry.Item3] = this.constrainedDofs[entry.Item1, entry.Item2];
            }
            elementToGlobalFreeDofs = globalFreeDofs;
            elementToGlobalConstrainedDofs = globalConstrainedDofs;
        }

        /// <summary>
        /// Index i = element local dof. Dictionary[i] = global dof.
        /// </summary>
        /// <param name="element"></param>
        /// <returns></returns>
        public IReadOnlyDictionary<int, int> MatchElementToGlobalEnrichedDofsOf(XContinuumElement2D element) //TODO: this can be an array. 
        {
            var globalEnrichedDofs = new Dictionary<int, int>();
            int elementDof = 0;
            foreach (XNode2D node in element.Nodes)
            {
                // TODO: Perhaps the nodal dof types should be decided by the element type (structural, continuum) 
                // and drawn from XXContinuumElement2D instead of from enrichment items
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    // TODO: Perhaps I could iterate directly on the dofs, ignoring dof types for performance, if the order is guarranteed
                    foreach (EnrichedDof dofType in enrichment.Dofs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        globalEnrichedDofs[elementDof++] = this.enrichedDofs[node, dofType];
                    }
                }
            }
            return globalEnrichedDofs;
        }

        public void ReorderUnconstrainedDofs(IReadOnlyList<int> permutation, bool oldToNew)
        {
            freeDofs.Reorder(permutation, oldToNew);
            enrichedDofs.Reorder(permutation, oldToNew);
        }

        public void WriteToConsole()
        {
            Console.WriteLine("Standard free dofs: node, dof, number");
            Console.WriteLine(freeDofs);
            Console.WriteLine("Enriched free dofs: node, dof, number");
            Console.WriteLine(enrichedDofs);
            Console.WriteLine("Standard constrained dofs: node, dof, number");
            Console.WriteLine(constrainedDofs);
        }
    }
}
