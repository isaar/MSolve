using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    class DOFEnumeratorInterleaved: DOFEnumeratorBase
    {
        private DOFEnumeratorInterleaved(int constrainedDofsCount, Table<XNode2D, DisplacementDOF, int> constrainedDofs,
            int enrichedDofsCount, Table<XNode2D, EnrichedDOF, int> enrichedDofs,
            int freeDofsCount, Table<XNode2D, DisplacementDOF, int> freeDofs) :
            base(constrainedDofsCount, constrainedDofs, enrichedDofsCount, enrichedDofs, freeDofsCount, freeDofs)
        {
        }

        public static DOFEnumeratorInterleaved Create(Model2D model)
        {
            // TODO: I should probably have a Constraint or Constraints class, to decouple this class from the collections 
            // used to represent constraints
            (Table<XNode2D, DisplacementDOF, int> freeDofs, Table<XNode2D, EnrichedDOF, int> enrichedDofs) = 
                EnumerateFreeEnrichedDofs(model);
            (int constrainedDofsCount, Table<XNode2D, DisplacementDOF, int> constrainedDofs) =
                EnumerateConstrainedDofs(model.Constraints);

            #region DEBUG code
            //Console.WriteLine("------------------------ DEBUG ------------------------------");
            //Console.WriteLine("Free standard dofs: ");
            //Console.WriteLine(freeDofs);
            //Console.WriteLine("Enriched dofs: ");
            //Console.WriteLine(enrichedDofs);
            //Console.WriteLine("Constrained dofs: ");
            //Console.WriteLine(constrainedDofs);
            //Console.WriteLine("------------------------ /DEBUG ------------------------------");
            #endregion

            return new DOFEnumeratorInterleaved(constrainedDofsCount, constrainedDofs, enrichedDofs.EntryCount, enrichedDofs,
                freeDofs.EntryCount, freeDofs);
        }

        private static (int constrainedDofsCount, Table<XNode2D, DisplacementDOF, int> constrainedDofs) EnumerateConstrainedDofs(
            ITable<XNode2D, DisplacementDOF, double> constraints)
        {
            var constrainedDofs = new Table<XNode2D, DisplacementDOF, int>();
            int dofCounter = 0;
            foreach (Tuple<XNode2D, DisplacementDOF, double> entry in constraints)
            {
                constrainedDofs[entry.Item1, entry.Item2] = dofCounter++;
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
        private static (Table<XNode2D, DisplacementDOF, int> freeDofs, Table<XNode2D, EnrichedDOF, int> enrichedDofs) 
            EnumerateFreeEnrichedDofs(Model2D model)
        {
            ITable<XNode2D, DisplacementDOF, double> constraints = model.Constraints;
            var freeDofs = new Table<XNode2D, DisplacementDOF, int>();
            var enrichedDofs = new Table<XNode2D, EnrichedDOF, int>();
            int dofCounter = 0;
            foreach (XNode2D node in model.Nodes)
            {
                // Standard free dofs. No rotational dofs. They can be X or Y. One or both of them may be constrained. 
                if (!constraints.Contains(node, DisplacementDOF.X)) freeDofs[node, DisplacementDOF.X] = dofCounter++;
                if (!constraints.Contains(node, DisplacementDOF.Y)) freeDofs[node, DisplacementDOF.Y] = dofCounter++;

                // Enriched dofs. No rotational dofs. They cannot be constrained.
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (EnrichedDOF dofType in enrichment.DOFs) enrichedDofs[node, dofType] = dofCounter++;
                }
            }
            return (freeDofs, enrichedDofs);
            //TODO: Also return each table's count, by keeping separate counters. This avoids the O(numRows) Table.Count()
        }
    }
}
