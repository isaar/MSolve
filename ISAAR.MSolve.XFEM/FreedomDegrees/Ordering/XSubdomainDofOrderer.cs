using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    //TODO: perhaps the subdomain should do all these itself
    class XSubdomainDofOrderer
    {
        protected readonly DofTable<EnrichedDof> enrichedDofs;

        private XSubdomainDofOrderer(int numEnrichedDofs, DofTable<EnrichedDof> enrichedDofs)
        {
            this.NumEnrichedDofs = numEnrichedDofs;
            this.enrichedDofs = enrichedDofs;
        }

        public int NumEnrichedDofs { get; }

        public static XSubdomainDofOrderer CreateNodeMajor(SortedSet<XNode2D> nodes)
        {
            var enrichedDofs = new DofTable<EnrichedDof>();
            int dofCounter = 0;
            foreach (XNode2D node in nodes)
            {
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (EnrichedDof dofType in enrichment.Dofs) enrichedDofs[node, dofType] = dofCounter++;
                }
            }
            return new XSubdomainDofOrderer(dofCounter, enrichedDofs);
        }

        public int GetDofIdx(XNode2D node, EnrichedDof dof)
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
                    foreach (EnrichedDof dofType in enrichment.Dofs)
                    {
                        subdomainEnrichedDofs[elementDof++] = this.enrichedDofs[node, dofType];
                    }
                }
            }
            return subdomainEnrichedDofs;
        }
    }
}
