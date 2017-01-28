using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Jump;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Enrichments
{
    class InterfaceEnrichment2D
    {
        private readonly ICurve2D curve;
        public IJumpEnrichment EnrichmentFunction { get; }
        public IReadOnlyDictionary<Node2D, double> EnrichedNodalDistances { get; }
        public IReadOnlyList<Element2D> IntersectedElements { get; }
        
        public InterfaceEnrichment2D(ICurve2D curve, IJumpEnrichment enrichmentFunction, 
            IEnumerable<Element2D> intersectedElements)
        {
            this.curve = curve;
            this.EnrichmentFunction = enrichmentFunction;

            var uniqueElements = new HashSet<Element2D>(intersectedElements);
            this.IntersectedElements = new List<Element2D>(uniqueElements);
            this.EnrichedNodalDistances = FindEnrichedNodalDistances();
        }

        private IReadOnlyDictionary<Node2D, double> FindEnrichedNodalDistances()
        {
            var uniqueNodes = new SortedSet<Node2D>();
            foreach (var element in IntersectedElements)
            {
                foreach (var node in element.Nodes)
                {
                    uniqueNodes.Add(node);
                }
            }

            var enrichedNodalDistances = new SortedDictionary<Node2D, double>();
            foreach (var node in uniqueNodes)
            {
                enrichedNodalDistances[node] = curve.SignedDistanceOf(node);
            }

            return enrichedNodalDistances;
        }
    }
}
