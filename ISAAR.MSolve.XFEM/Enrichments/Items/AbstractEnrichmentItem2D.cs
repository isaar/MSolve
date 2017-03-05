using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class AbstractEnrichmentItem2D: IEnrichmentItem2D
    {
        protected List<XElement2D> affectedElements;

        public int ID { get; }
        public IReadOnlyList<IEnrichmentFunction2D> EnrichmentFunctions { get; protected set; }
        public IReadOnlyList<ArtificialDOFType> DOFs { get; protected set; }
        public IReadOnlyList<XElement2D> AffectedElements { get { return affectedElements; } }

        protected AbstractEnrichmentItem2D()
        {
            this.affectedElements = new List<XElement2D>();
        }

        public void AffectElement(XElement2D element)
        {
            // TODO: There should be a check here or this method should be private.
            if (!affectedElements.Contains(element)) affectedElements.Add(element);
        }

        //ERROR: This method overwrites the enrichment functions from previous items. Thus if there is both a tip and a heaviside, only the last will be saved
        //TODO: This method should also handle the removal of enrichments form a node.
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
                node.EnrichmentItems = new IEnrichmentItem2D[] { this };
                node.EnrichmentFunctions = allEnrichments;
            }
        }
    }
}
