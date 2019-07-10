using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class MaterialInterface2D : AbstractEnrichmentItem2D
    {
        public enum Subdomain { Positive, Negative, Boundary }   
             
        private readonly RampFunction2D enrichmentFunction;

        public MaterialInterface2D(ICurve2D geometry)
        {
            this.Discontinuity = geometry;
            this.enrichmentFunction = new RampFunction2D();
            this.Dofs = new EnrichedDof[] {
                new EnrichedDof(enrichmentFunction, StructuralDof.TranslationX),
                new EnrichedDof(enrichmentFunction, StructuralDof.TranslationY)
            };
            this.ElementIntersections = new Dictionary<XContinuumElement2D, CartesianPoint[]>();
        }

        public ICurve2D Discontinuity { get; }
        public Dictionary<XContinuumElement2D, CartesianPoint[]> ElementIntersections { get; }

        public override double[] EvaluateFunctionsAt(XNode node)
        {
            double signedDistance = Discontinuity.SignedDistanceOf(node);
            return new double[] { enrichmentFunction.EvaluateAt(signedDistance) };
        }

        public override EvaluatedFunction2D[] EvaluateAllAt(NaturalPoint point, XContinuumElement2D element,
             EvalInterpolation2D interpolation)
        {
            CartesianPoint cartesianPoint = interpolation.TransformPointNaturalToGlobalCartesian();
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
        public Subdomain LocatePoint(CartesianPoint point)
        {
            int sign = Math.Sign(Discontinuity.SignedDistanceOf(point));
            if (sign < 0) return Subdomain.Negative;
            else if (sign > 0) return Subdomain.Positive;
            else return Subdomain.Boundary;
        }

        public override IReadOnlyList<CartesianPoint> IntersectionPointsForIntegration(XContinuumElement2D element)
        {
            CartesianPoint[] intersectionPoints;
            bool alreadyIntersected = ElementIntersections.TryGetValue(element, out intersectionPoints);
            if (alreadyIntersected) return intersectionPoints;
            else
            {
                throw new NotImplementedException();
            }
        }
    }
}
