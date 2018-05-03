using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    //TODO: perhaps the subdomain should do all these itself
    class DOFEnumeratorXSubdomain
    {
        protected readonly DOFTable<EnrichedDOF> enrichedDofs;

        private DOFEnumeratorXSubdomain(int numEnrichedDofs, DOFTable<EnrichedDOF> enrichedDofs)
        {
            this.NumEnrichedDofs = numEnrichedDofs;
            this.enrichedDofs = enrichedDofs;
        }

        public int NumEnrichedDofs { get; }

        public static DOFEnumeratorXSubdomain CreateNodeMajor(SortedSet<XNode2D> nodes)
        {
            var enrichedDofs = new DOFTable<EnrichedDOF>();
            int dofCounter = 0;
            foreach (XNode2D node in nodes)
            {
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (EnrichedDOF dofType in enrichment.DOFs) enrichedDofs[node, dofType] = dofCounter++;
                }
            }
            return new DOFEnumeratorXSubdomain(dofCounter, enrichedDofs);
        }

        public int GetDofIdx(XNode2D node, EnrichedDOF dof)
        {
            return enrichedDofs[node, dof];
        }

        public IEnumerable<int> GetDofIndices(XNode2D node)
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

        //TODO: this could be an array
        public IReadOnlyDictionary<int, int> MatchElementToSubdomainEnrichedDofs(XContinuumElement2D element)
        {
            var subdomainEnrichedDofs = new Dictionary<int, int>();
            int elementDof = 0;
            foreach (XNode2D node in element.Nodes)
            {
                // TODO: Perhaps the nodal dof types should be decided by the element type (structural, continuum) 
                // and drawn from XXContinuumElement2D instead of from enrichment items
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    // TODO: Perhaps I could iterate directly on the dofs, ignoring dof types for performance, if the order is guarranteed
                    foreach (EnrichedDOF dofType in enrichment.DOFs)
                    {
                        subdomainEnrichedDofs[elementDof++] = this.enrichedDofs[node, dofType];
                    }
                }
            }
            return subdomainEnrichedDofs;
        }
    }
}
