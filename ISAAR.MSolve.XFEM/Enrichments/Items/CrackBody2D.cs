using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class CrackBody2D : IEnrichmentItem2D
    {
        private List<XElement2D> splitElements;

        public ICurve2D Discontinuity { get; }
        public IReadOnlyList<IEnrichmentFunction2D> EnrichmentFunctions { get; }
        public IReadOnlyList<XElement2D> AffectedElements { get { return splitElements; } }

        //TODO: allow the user to specify which Heaviside function to pass. 
        // First specialize IEnrichmentFunction2D into IDiscontinuityFunction2D to allow only Heaviside, Sign or 
        // Heaviside approximation functions.
        public CrackBody2D(ICurve2D discontinuity) 
        {
            this.splitElements = new List<XElement2D>();
            this.Discontinuity = discontinuity;
            this.EnrichmentFunctions = new IEnrichmentFunction2D[] { new SignFunction2D(this) };
        }

        public void AffectElement(XElement2D element)
        {
            // TODO: There should be a check here or this method should be private.
            if (!splitElements.Contains(element)) splitElements.Add(element);
        }

        // Ok for the 1st time. What about updates when only some enrichments must be cleared/changed?
        public void EnrichNodes()
        {
            // Find all unique affected nodes.
            HashSet<XNode2D> nodes = new HashSet<XNode2D>();
            foreach (var element in AffectedElements) nodes.UnionWith(element.Nodes);
            
            foreach (var node in nodes)
            {
                var allEnrichments = new Tuple<IEnrichmentFunction2D, double>[EnrichmentFunctions.Count];
                int enrichmentCounter = 0;
                foreach (var enrichmentFunction in EnrichmentFunctions)
                {
                    allEnrichments[enrichmentCounter] = 
                        new Tuple<IEnrichmentFunction2D, double>(enrichmentFunction, enrichmentFunction.EvalueAt(node));
                    ++enrichmentCounter;
                }
                node.EnrichmentFunctions = allEnrichments;
            }
        }
    }
}
