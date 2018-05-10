using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.Descriptions;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class MaterialInterface2D : AbstractEnrichmentItem2D
    {
        public enum Subdomain { Positive, Negative, Boundary }   
             
        public IGeometryDescription2D Discontinuity { get; }
        private readonly RampFunction2D enrichmentFunction;

        public MaterialInterface2D(IGeometryDescription2D geometry)
        {
            this.Discontinuity = geometry;
            this.enrichmentFunction = new RampFunction2D();
            this.Dofs = new EnrichedDof[] {
                new EnrichedDof(enrichmentFunction, DisplacementDof.X),
                new EnrichedDof(enrichmentFunction, DisplacementDof.Y)
            };
        }

        public override double[] EvaluateFunctionsAt(XNode2D node)
        {
            double signedDistance = Discontinuity.SignedDistanceOf(node);
            return new double[] { enrichmentFunction.EvaluateAt(signedDistance) };
        }

        public override EvaluatedFunction2D[] EvaluateAllAt(INaturalPoint2D point, XContinuumElement2D element,
             EvaluatedInterpolation2D interpolation)
        {
            ICartesianPoint2D cartesianPoint = interpolation.TransformPointNaturalToGlobalCartesian(point);
            double signedDistance = Discontinuity.SignedDistanceOf(cartesianPoint);
            Vector2 normalVector = Discontinuity.NormalVectorThrough(cartesianPoint);
            return new EvaluatedFunction2D[] { enrichmentFunction.EvaluateAllAt(signedDistance, normalVector) };
        }

        // TODO: add some tolerance when checking around 0. Perhaps all this is not needed though and I could even 
        // ignore points on the interface. It certainly needs a better name
        /// <summary>
        /// Finds the subdomain where the requested cartesian point lies.
        /// </summary>
        /// <param name="point"></param>
        /// <param name="subdomain">The posi</param>
        /// <returns></returns>
        public Subdomain LocatePoint(ICartesianPoint2D point)
        {
            int sign = Math.Sign(Discontinuity.SignedDistanceOf(point));
            if (sign < 0) return Subdomain.Negative;
            else if (sign > 0) return Subdomain.Positive;
            else return Subdomain.Boundary;
        }

        public override IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element)
        {
            return Discontinuity.IntersectionWith(element);
        }
    }
}
