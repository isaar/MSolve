using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.Descriptions;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class CrackBody2D : AbstractEnrichmentItem2D
    {
        public IGeometryDescription2D Discontinuity { get; }
        private readonly SignFunction2D enrichmentFunction;

        //TODO: allow the user to specify which Heaviside function to pass. 
        // First specialize IEnrichmentFunction2D into IDiscontinuityFunction2D to allow only Heaviside, Sign or 
        // Heaviside approximation functions.
        public CrackBody2D(IGeometryDescription2D discontinuity)
        {
            this.Discontinuity = discontinuity;
            enrichmentFunction = new SignFunction2D();
            this.DOFs = new ArtificialDOFType[] {
                new ArtificialDOFType(enrichmentFunction, StandardDOFType.X),
                new ArtificialDOFType(enrichmentFunction, StandardDOFType.Y)
            };
        }

        public override double[] EvaluateFunctionsAt(ICartesianPoint2D point)
        {
            double signedDistance = Discontinuity.SignedDistanceOf(point);
            return new double[] { enrichmentFunction.EvaluateAt(signedDistance) };
        }

        public override EvaluatedFunction2D[] EvaluateAllAt(INaturalPoint2D point, IReadOnlyList<XNode2D> elementNodes,
             EvaluatedInterpolation2D interpolation)
        {
            ICartesianPoint2D cartesianPoint = interpolation.TransformPointNaturalToGlobalCartesian(point);
            double signedDistance = Discontinuity.SignedDistanceOf(cartesianPoint);
            return new EvaluatedFunction2D[] { enrichmentFunction.EvaluateAllAt(signedDistance) };
        }

        public override IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element)
        {
            var uniquePoints = new HashSet<ICartesianPoint2D>(element.Nodes);
            foreach (ICartesianPoint2D point in Discontinuity.IntersectionWith(element))
            {
                uniquePoints.Add(point);
            }
            return new List<ICartesianPoint2D>(uniquePoints);
        }
    }
}
