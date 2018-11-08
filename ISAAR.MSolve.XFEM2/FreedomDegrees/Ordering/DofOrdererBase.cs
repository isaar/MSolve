using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

//TODO: I should use interchangeble builders, rather than subclasses.
//TODO: constrained dofs should be ignored (set to -1) and the effect of constraints should be done via equivalent forces.
namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    class DofOrdererBase: IDofOrderer
    {
        protected readonly DofTable<DisplacementDof> constrainedDofs;
        protected readonly DofTable<EnrichedDof> enrichedDofs;
        protected readonly DofTable<DisplacementDof> standardDofs;

        public int NumConstrainedDofs { get; }
        public int NumEnrichedDofs { get; }
        public int NumStandardDofs { get; }

        protected DofOrdererBase(int constrainedDofsCount, DofTable<DisplacementDof> constrainedDofs, 
            int enrichedDofsCount, DofTable<EnrichedDof> enrichedDofs,
            int standardDofsCount, DofTable<DisplacementDof> standardDofs)
        {
            this.NumConstrainedDofs = constrainedDofsCount;
            this.constrainedDofs = constrainedDofs;
            this.NumEnrichedDofs = enrichedDofsCount;
            this.enrichedDofs = enrichedDofs;
            this.NumStandardDofs = standardDofsCount;
            this.standardDofs = standardDofs;
        }

        public IDofOrderer DeepCopy()
        {
            return new DofOrdererBase(NumConstrainedDofs, constrainedDofs.DeepCopy(),
                NumEnrichedDofs, enrichedDofs.DeepCopy(), NumStandardDofs, standardDofs.DeepCopy());
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
            foreach ((XNode2D node, DisplacementDof dofType, int dofIdx) in elementDofs)
            {
                bool isStandard = this.standardDofs.TryGetValue(node, dofType, out int globalStandardDof);
                if (isStandard) elementVector[dofIdx] = globalFreeVector[globalStandardDof];
                else
                {
                    int globalConstrainedDof = this.constrainedDofs[node, dofType];
                    elementVector[dofIdx] = globalConstrainedVector[globalConstrainedDof];
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
            DofTable<EnrichedDof> elementDofs = element.GetEnrichedDofs();
            double[] elementVector = new double[elementDofs.EntryCount];
            foreach ((XNode2D node, EnrichedDof dofType, int dofIdx) in elementDofs)
            {
                int globalEnrichedDof = enrichedDofs[node, dofType];
                elementVector[dofIdx] = globalFreeVector[globalEnrichedDof];
            }
            return Vector.CreateFromArray(elementVector);
        }

        public ITable<XNode2D, EnrichedDof, double> GatherEnrichedNodalDisplacements(Model2D model, Vector solution)
        {
            var table = new Table<XNode2D, EnrichedDof, double>();
            foreach (var row in enrichedDofs)
            {
                table[row.row, row.col] = solution[row.val];
            }
            return table;
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

        public List<int> GetConstrainedDofsOf(XContinuumElement2D element)
        {
            var globalDofs = new List<int>();
            foreach (var nodeDofLocal in element.GetStandardDofs())
            {
                bool isConstrained = constrainedDofs.TryGetValue(nodeDofLocal.row, nodeDofLocal.col, out int globalDof);
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
                globalDofs.Add(enrichedDofs[nodeDofLocal.row, nodeDofLocal.col]); // It must be included.
            }
            return globalDofs;
        }

        // Would it be faster to return the Dictionary<StandardDofType, int> for consecutive accesses of the dofs of this node? 
        // Dictionary's Read operation is supposed to be O(1), but still...
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
                bool isStandard = standardDofs.TryGetValue(nodeDofLocal.row, nodeDofLocal.col, out int globalDof);
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
            foreach ((XNode2D node, DisplacementDof dofType, int dofIdx) in elementDofs)
            {
                bool isStandard = this.standardDofs.TryGetValue(node, dofType, out int standardGlobalDof);
                if (isStandard) globalStandardDofs[dofIdx] = standardGlobalDof;
                else globalConstrainedDofs[dofIdx] = this.constrainedDofs[node, dofType];
            }
            elementToGlobalStandardDofs = globalStandardDofs;
            elementToGlobalConstrainedDofs = globalConstrainedDofs;
        }

        public void ReorderUnconstrainedDofs(IReadOnlyList<int> permutation, bool oldToNew)
        {
            standardDofs.Reorder(permutation, oldToNew);
            enrichedDofs.Reorder(permutation, oldToNew);
        }

        public void WriteToConsole()
        {
            Console.WriteLine("Standard free dofs: node, dof, number");
            Console.WriteLine(standardDofs);
            Console.WriteLine("Enriched free dofs: node, dof, number");
            Console.WriteLine(enrichedDofs);
            Console.WriteLine("Standard constrained dofs: node, dof, number");
            Console.WriteLine(constrainedDofs);
        }
    }
}
