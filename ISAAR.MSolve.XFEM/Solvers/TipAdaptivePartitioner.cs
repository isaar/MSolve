using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

//TODO: At the start of each iteration, find which region the tip belongs to, instead of taking the current subdomain of the 
//      tip element. Then adapt the elements/nodes/subdomains around the crack tip.
//TODO: It is possible that due to this exchange of elements between subdomains, some corner nodes change. How to handle this?
namespace ISAAR.MSolve.XFEM.Solvers
{
    public class TipAdaptivePartitioner
    {
        private readonly ICrackDescription crack;
        //private readonly IReadOnlyList<IRegion2D> initialRegions; //TODO: this should be used to make sure the crack does not end up growing its initial subdomain too much.

        public TipAdaptivePartitioner(ICrackDescription crack)
        {
            this.crack = crack;
        }

        public HashSet<ISubdomain> UpdateSubdomains()
        {
            var modifiedSubdomains = new HashSet<XSubdomain>();

            //TODO: this can be done without accessing the single crack. We just need to access the same tip element and nodes
            foreach (ISingleCrack singleCrack in crack.SingleCracks)
            {
                // Find the subdomain that contains the crack tip
                CartesianPoint tip = singleCrack.CrackTips[0];
                IXFiniteElement tipElement = singleCrack.CrackTipElements[tip][0]; // If there are more neighboring ones, we don't need them
                XSubdomain tipSubdomain = tipElement.Subdomain;

                // Find all elements that have at least one tip enriched node
                //TODO: this can be merged with the next subtask
                var tipEnrichedElements = new HashSet<IXFiniteElement>();
                ISet<XNode> tipNodes = singleCrack.CrackTipNodesNew[singleCrack.CrackTipEnrichments];
                foreach (XNode node in tipNodes) tipEnrichedElements.UnionWith(node.ElementsDictionary.Values);

                // Find which of these elements must be removed and from which subdomains
                foreach (IXFiniteElement element in tipEnrichedElements)
                {
                    if (element.Subdomain != tipSubdomain)
                    {
                        //TODO: For once, having a Dictionary would be more performant. Still this could be implemented as as 
                        //      binary search.
                        element.Subdomain.Elements.Remove(element); 
                        modifiedSubdomains.Add(element.Subdomain);

                        //TODO: It should be added to a position that will not affect the bandwidth to adversely.
                        tipSubdomain.Elements.Add(element); 
                        modifiedSubdomains.Add(tipSubdomain);
                        element.Subdomain = tipSubdomain;
                    }
                }

                // Update the subdomains of each modified node
                var modifiedNodes = new HashSet<XNode>();
                foreach (IXFiniteElement element in tipEnrichedElements) modifiedNodes.UnionWith(element.Nodes);
                foreach (XNode node in modifiedNodes)
                {
                    node.SubdomainsDictionary.Clear();
                    node.BuildXSubdomainDictionary();
                }
            }

            // Find the new nodes in the subdomains with added or removed elements
            foreach (XSubdomain subdomain in modifiedSubdomains) subdomain.DefineNodesFromElements();

            return new HashSet<ISubdomain>(modifiedSubdomains);
        }
    }
}
