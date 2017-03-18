using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    abstract class AbstractEnrichmentItem2D: IEnrichmentItem2D
    {
        protected List<XContinuumElement2D> affectedElements;

        public int ID { get; }
        public IReadOnlyList<IEnrichmentFunction2D> EnrichmentFunctions { get; protected set; }
        public IReadOnlyList<ArtificialDOFType> DOFs { get; protected set; }
        public IReadOnlyList<XContinuumElement2D> AffectedElements { get { return affectedElements; } }

        protected AbstractEnrichmentItem2D()
        {
            this.affectedElements = new List<XContinuumElement2D>();
        }

        public void AffectElement(XContinuumElement2D element)
        {
            // TODO: There should be a check here or this method should be private.
            if (!affectedElements.Contains(element))
            {
                affectedElements.Add(element);
                element.EnrichmentItems.Add(this);
            }
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
                node.EnrichmentItems = new List<IEnrichmentItem2D> { };
                node.EnrichmentItems.Add(this);

                var enrichments = new List<Tuple<IEnrichmentFunction2D, double>>(EnrichmentFunctions.Count);
                foreach (var function in EnrichmentFunctions)
                {
                    enrichments.Add(new Tuple<IEnrichmentFunction2D, double>(function, function.EvalueAt(node)));
                }
                node.EnrichmentFunctions = enrichments;
            }
        }

        public void EnrichNode(XNode2D node) // TODO: delete after debugging
        {
            node.EnrichmentItems.Add(this);
            foreach (var function in EnrichmentFunctions)
            {
                node.EnrichmentFunctions.Add(
                    new Tuple<IEnrichmentFunction2D, double>(function, function.EvalueAt(node)));
            }
        }

        public void EnrichElement(XContinuumElement2D element)
        {
            element.EnrichmentItems.Add(this);
        }

        public abstract IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element);
    }
}
