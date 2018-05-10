using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
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

        public void EnrichNode(XNode2D node) // TODO: this should not be done manually
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

        public abstract double[] EvaluateFunctionsAt(XNode2D node);
        public abstract EvaluatedFunction2D[] EvaluateAllAt(INaturalPoint2D point, XContinuumElement2D element,
             EvaluatedInterpolation2D interpolation);

        public abstract IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element);
    }
}
