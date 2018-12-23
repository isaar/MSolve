using System;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    class InterleavedDofOrderer: DofOrdererBase
    {
        private InterleavedDofOrderer(int constrainedDofsCount, DofTable<DisplacementDof> constrainedDofs,
            int enrichedDofsCount, DofTable<EnrichedDof> enrichedDofs,
            int standardDofsCount, DofTable<DisplacementDof> standardDofs) :
            base(constrainedDofsCount, constrainedDofs, enrichedDofsCount, enrichedDofs, standardDofsCount, standardDofs)
        {
        }

        public static InterleavedDofOrderer Create(Model2D model)
        {
            // TODO: I should probably have a Constraint or Constraints class, to decouple this class from the collections 
            // used to represent constraints
            (DofTable<DisplacementDof> standardDofs, DofTable<EnrichedDof> enrichedDofs) = 
                OrderUnconstrainedDofs(model);
            (int constrainedDofsCount, DofTable<DisplacementDof> constrainedDofs) =
                OrderConstrainedDofs(model.Constraints);

            #region DEBUG code
            //Console.WriteLine("------------------------ DEBUG ------------------------------");
            //Console.WriteLine("Free standard dofs: ");
            //Console.WriteLine(standardDofs);
            //Console.WriteLine("Enriched dofs: ");
            //Console.WriteLine(enrichedDofs);
            //Console.WriteLine("Constrained dofs: ");
            //Console.WriteLine(constrainedDofs);
            //Console.WriteLine("------------------------ /DEBUG ------------------------------");
            #endregion

            return new InterleavedDofOrderer(constrainedDofsCount, constrainedDofs, enrichedDofs.EntryCount, enrichedDofs,
                standardDofs.EntryCount, standardDofs);
        }

        private static (int constrainedDofsCount, DofTable<DisplacementDof> constrainedDofs) OrderConstrainedDofs(
            ITable<XNode2D, DisplacementDof, double> constraints)
        {
            var constrainedDofs = new DofTable<DisplacementDof>();
            int dofCounter = 0;
            foreach ((XNode2D node, DisplacementDof dofType, double displacement) in constraints)
            {
                constrainedDofs[node, dofType] = dofCounter++;
            }
            return (dofCounter, constrainedDofs);
        }

        /// <summary>
        /// Ordering is node major, then enrichment major and finally axis minor: [... uNx, uNy, aN1x, aN1y, aN2x, aN2y, ...], 
        /// where N is the node index. Warning: This only works for continuum mechanics, where there are no rotational dofs at 
        /// certain nodes.
        /// </summary>
        /// <param name="model"></param>
        /// <returns></returns>
        private static (DofTable<DisplacementDof> standardDofs, DofTable<EnrichedDof> enrichedDofs) 
            OrderUnconstrainedDofs(Model2D model)
        {
            ITable<XNode2D, DisplacementDof, double> constraints = model.Constraints;
            var standardDofs = new DofTable<DisplacementDof>();
            var enrichedDofs = new DofTable<EnrichedDof>();
            int dofCounter = 0;
            foreach (XNode2D node in model.Nodes)
            {
                // Standard free dofs. No rotational dofs. They can be X or Y. One or both of them may be constrained. 
                if (!constraints.Contains(node, DisplacementDof.X)) standardDofs[node, DisplacementDof.X] = dofCounter++;
                if (!constraints.Contains(node, DisplacementDof.Y)) standardDofs[node, DisplacementDof.Y] = dofCounter++;

                // Enriched dofs. No rotational dofs. They cannot be constrained.
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (EnrichedDof dofType in enrichment.Dofs) enrichedDofs[node, dofType] = dofCounter++;
                }
            }
            return (standardDofs, enrichedDofs);
            //TODO: Also return each table's count, by keeping separate counters. This avoids the O(numRows) Table.Count()
        }
    }
}
