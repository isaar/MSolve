using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    //TODO: perhaps the subdomain should do all these itself
    class XSubdomainDofOrderer
    {
        //TODO: perhaps the rows, columns of the 2 tables can be shared, for less memory and faster simultaneous iteration
        protected readonly DofTable<EnrichedDof> globalEnrichedDofs; 
        protected readonly DofTable<EnrichedDof> subdomainEnrichedDofs;

        private XSubdomainDofOrderer(int numEnrichedDofs, DofTable<EnrichedDof> subdomainEnrichedDofs,
            DofTable<EnrichedDof> globalEnrichedDofs)
        {
            this.NumEnrichedDofs = numEnrichedDofs;
            this.subdomainEnrichedDofs = subdomainEnrichedDofs;
            this.globalEnrichedDofs = globalEnrichedDofs;
        }

        public int NumEnrichedDofs { get; }

        public static XSubdomainDofOrderer CreateNodeMajor(XSubdomain2D subdomain, int globalIndicesStart) //TODO: also add AMD reordering
        {
            var subdomainEnrichedDofs = new DofTable<EnrichedDof>();
            var globalEnrichedDofs = new DofTable<EnrichedDof>();
            int dofCounter = 0;
            foreach (XNode2D node in subdomain.AllNodes)
            {
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (EnrichedDof dofType in enrichment.Dofs)
                    {
                        subdomainEnrichedDofs[node, dofType] = dofCounter;
                        globalEnrichedDofs[node, dofType] = globalIndicesStart + dofCounter; // Then I could just subtract the offest to go local -> global
                        ++dofCounter;
                    }
                }
            }
            return new XSubdomainDofOrderer(dofCounter, subdomainEnrichedDofs, globalEnrichedDofs);
        }

        //public int GetDofIdx(XNode2D node, EnrichedDof dof)
        //{
        //    return subdomainEnrichedDofs[node, dof];
        //}

        //public IEnumerable<int> GetEnrichedDofs(XNode2D node)
        //{
        //    try
        //    {
        //        return subdomainEnrichedDofs.GetValuesOfRow(node);
        //    }
        //    catch (KeyNotFoundException)
        //    {
        //        return new int[] { }; // TODO: It would be better if this method is not called for a non enriched node.
        //    }
        //}

        public (IReadOnlyDictionary<int, int> element2Subdomain, IReadOnlyDictionary<int, int> element2Global) 
            MatchElementToSubdomainAndGlobalEnrichedDofs(XContinuumElement2D element)
        {
            var element2Subdomain = new Dictionary<int, int>();
            var element2Global = new Dictionary<int, int>();
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
                        element2Subdomain[elementDof] = this.subdomainEnrichedDofs[node, dofType];
                        element2Global[elementDof] = this.globalEnrichedDofs[node, dofType];
                        ++elementDof;
                    }
                }
            }
            return (element2Subdomain, element2Global);
        }

        public void ReorderEnrichedSubdomainDofs(OrderingAMD orderingAlgorithm)
        {
            throw new NotImplementedException();
        }

        public void WriteToConsole()
        {
            Console.WriteLine("Enriched dofs - subdomain order: ");
            Console.WriteLine(subdomainEnrichedDofs);
            Console.WriteLine("Enriched dofs - global order: ");
            Console.WriteLine(globalEnrichedDofs);
        }
    }
}
