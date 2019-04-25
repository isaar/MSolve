using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Utilities;

// TODO: this class should not be associated with the whole crack geometry, just the part that stores a single branch.
namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    public class CrackBodyEnrichment2D : IEnrichmentItem2D
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
                new EnrichedDof(enrichmentFunction, StructuralDof.TranslationX),
                new EnrichedDof(enrichmentFunction, StructuralDof.TranslationY)
            };
        }
        
        public double[] EvaluateFunctionsAt(XNode node)
        {
            double signedDistance = crackDescription.SignedDistanceOf(node);
            return new double[] { enrichmentFunction.EvaluateAt(signedDistance) };
        }

        public EvaluatedFunction2D[] EvaluateAllAt(NaturalPoint point, XContinuumElement2D element,
             EvalInterpolation2D interpolation)
        {
            CartesianPoint cartesianPoint = interpolation.TransformPointNaturalToGlobalCartesian();
            double signedDistance = crackDescription.SignedDistanceOf(point, element, interpolation);
            return new EvaluatedFunction2D[] { enrichmentFunction.EvaluateAllAt(signedDistance) };
        }


        // TODO: delete this
        public IReadOnlyList<CartesianPoint> IntersectionPointsForIntegration(XContinuumElement2D element)
        {
            throw new NotImplementedException();
        }
    }
}
