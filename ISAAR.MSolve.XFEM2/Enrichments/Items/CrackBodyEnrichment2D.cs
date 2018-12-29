using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;

// TODO: this class should not be associated with the whole crack geometry, just the part that stores a single branch.
namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class CrackBodyEnrichment2D : IEnrichmentItem2D
    {
        private readonly IHeavisideFunction2D enrichmentFunction;
        public IReadOnlyList<EnrichedDof> Dofs { get; }
        public ISingleCrack crackDescription;

        public CrackBodyEnrichment2D(ISingleCrack crackDescription): this(crackDescription, new SignFunction2D())
        {
        }

        public CrackBodyEnrichment2D(ISingleCrack crackDescription, IHeavisideFunction2D enrichmentFunction)
        {
            this.crackDescription = crackDescription;
            this.enrichmentFunction = enrichmentFunction;
            this.Dofs = new EnrichedDof[] {
                new EnrichedDof(enrichmentFunction, DisplacementDof.X),
                new EnrichedDof(enrichmentFunction, DisplacementDof.Y)
            };
        }
        
        public double[] EvaluateFunctionsAt(XNode2D node)
        {
            double signedDistance = crackDescription.SignedDistanceOf(node);
            return new double[] { enrichmentFunction.EvaluateAt(signedDistance) };
        }

        public EvaluatedFunction2D[] EvaluateAllAt(INaturalPoint2D point, XContinuumElement2D element,
             EvaluatedInterpolation2D interpolation)
        {
            ICartesianPoint2D cartesianPoint = interpolation.TransformPointNaturalToGlobalCartesian(point);
            double signedDistance = crackDescription.SignedDistanceOf(point, element, interpolation);
            return new EvaluatedFunction2D[] { enrichmentFunction.EvaluateAllAt(signedDistance) };
        }


        // TODO: delete this
        public IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element)
        {
            throw new NotImplementedException();
        }
    }
}
