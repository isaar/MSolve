using System.Collections.Generic;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    abstract class AbstractEnrichmentItem2D: IEnrichmentItem2D
    {
        protected List<XContinuumElement2D> affectedElements;

        public IReadOnlyList<EnrichedDof> Dofs { get; protected set; }
        public IReadOnlyList<XContinuumElement2D> AffectedElements { get { return affectedElements; } }

        protected AbstractEnrichmentItem2D()
        {
            this.affectedElements = new List<XContinuumElement2D>();
        }

        public void EnrichNode(XNode node) // TODO: this should not be done manually
        {
            double[] enrichmentValues = EvaluateFunctionsAt(node);
            node.EnrichmentItems.Add(this, enrichmentValues);
        }

        public void EnrichElement(XContinuumElement2D element)
        {
            if (!affectedElements.Contains(element))
            {
                affectedElements.Add(element);
                element.EnrichmentItems.Add(this);
            }
        }

        public abstract double[] EvaluateFunctionsAt(XNode node);
        public abstract EvaluatedFunction2D[] EvaluateAllAt(NaturalPoint point, XContinuumElement2D element,
             EvalInterpolation2D interpolation);

        public abstract IReadOnlyList<CartesianPoint> IntersectionPointsForIntegration(XContinuumElement2D element);
    }
}
