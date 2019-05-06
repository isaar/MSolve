using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    class InterleavedDofOrderer: DofOrdererBase
    {
        private InterleavedDofOrderer(int constrainedDofsCount, DofTable<StructuralDof> constrainedDofs,
            int enrichedDofsCount, DofTable<EnrichedDof> enrichedDofs,
            int standardDofsCount, DofTable<StructuralDof> standardDofs) :
            base(constrainedDofsCount, constrainedDofs, enrichedDofsCount, enrichedDofs, standardDofsCount, standardDofs)
        {
        }

        public static InterleavedDofOrderer Create(Model2D_old model)
        {
            // TODO: I should probably have a Constraint or Constraints class, to decouple this class from the collections 
            // used to represent constraints
            (DofTable<StructuralDof> standardDofs, DofTable<EnrichedDof> enrichedDofs) = 
                OrderUnconstrainedDofs(model);
            (int constrainedDofsCount, DofTable<StructuralDof> constrainedDofs) =
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

        private static (int constrainedDofsCount, DofTable<StructuralDof> constrainedDofs) OrderConstrainedDofs(
            ITable<XNode, StructuralDof, double> constraints)
        {
            var constrainedDofs = new DofTable<StructuralDof>();
            int dofCounter = 0;
            foreach ((XNode node, StructuralDof dofType, double displacement) in constraints)
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
        private static (DofTable<StructuralDof> standardDofs, DofTable<EnrichedDof> enrichedDofs) 
            OrderUnconstrainedDofs(Model2D_old model)
        {
            ITable<XNode, StructuralDof, double> constraints = model.Constraints;
            var standardDofs = new DofTable<StructuralDof>();
            var enrichedDofs = new DofTable<EnrichedDof>();
            int dofCounter = 0;
            foreach (XNode node in model.Nodes)
            {
                // Standard free dofs. No rotational dofs. They can be X or Y. One or both of them may be constrained. 
                if (!constraints.Contains(node, StructuralDof.TranslationX)) standardDofs[node, StructuralDof.TranslationX] = dofCounter++;
                if (!constraints.Contains(node, StructuralDof.TranslationY)) standardDofs[node, StructuralDof.TranslationY] = dofCounter++;

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
