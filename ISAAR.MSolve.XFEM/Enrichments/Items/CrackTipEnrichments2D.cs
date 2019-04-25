using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Utilities;

// TODO: this class should not be associated with the whole crack geometry, just the part that stores the crack tip.
// Also it should be immutable, rather than having a settable TipSystem property. 
namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    public class CrackTipEnrichments2D : IEnrichmentItem2D
    {
        private readonly IReadOnlyList<ITipFunction> enrichmentFunctions;
        private readonly ISingleCrack crackDescription;
        private readonly CrackTipPosition tipPosition;
        public IReadOnlyList<EnrichedDof> Dofs { get; }

        public CrackTipEnrichments2D(ISingleCrack crackDescription, CrackTipPosition tipPosition) : 
            this(crackDescription, tipPosition, new ITipFunction[] 
            {
                new IsotropicBrittleTipFunctions2D.Func1(),
                new IsotropicBrittleTipFunctions2D.Func2(),
                new IsotropicBrittleTipFunctions2D.Func3(),
                new IsotropicBrittleTipFunctions2D.Func4()
            })
        {
        }

        public CrackTipEnrichments2D(ISingleCrack crackDescription, CrackTipPosition tipPosition,
            IReadOnlyList<ITipFunction> enrichmentFunctions)
        {
            this.crackDescription = crackDescription;
            this.tipPosition = tipPosition;
            this.enrichmentFunctions = enrichmentFunctions;
            this.Dofs = new EnrichedDof[]
            {
                new EnrichedDof(enrichmentFunctions[0], StructuralDof.TranslationX),
                new EnrichedDof(enrichmentFunctions[0], StructuralDof.TranslationY),
                new EnrichedDof(enrichmentFunctions[1], StructuralDof.TranslationX),
                new EnrichedDof(enrichmentFunctions[1], StructuralDof.TranslationY),
                new EnrichedDof(enrichmentFunctions[2], StructuralDof.TranslationX),
                new EnrichedDof(enrichmentFunctions[2], StructuralDof.TranslationY),
                new EnrichedDof(enrichmentFunctions[3], StructuralDof.TranslationX),
                new EnrichedDof(enrichmentFunctions[3], StructuralDof.TranslationY),
            };

        }

        public TipCoordinateSystem TipSystem { get; set; }

        public double[] EvaluateFunctionsAt(XNode node)
        {
            PolarPoint2D polarPoint = TipSystem.TransformPointGlobalCartesianToLocalPolar(node);
            var enrichments = new double[enrichmentFunctions.Count];
            for (int i = 0; i < enrichments.Length; ++i)
            {
                enrichments[i] = enrichmentFunctions[i].EvaluateAt(polarPoint);
            }
            return enrichments;
        }

        public EvaluatedFunction2D[] EvaluateAllAt(NaturalPoint point, XContinuumElement2D element,
             EvalInterpolation2D interpolation)
        {
            PolarPoint2D polarPoint = TipSystem.TransformPointGlobalCartesianToLocalPolar(
                interpolation.TransformPointNaturalToGlobalCartesian());
            TipJacobians tipJacobians = TipSystem.CalculateJacobiansAt(polarPoint);

            var enrichments = new EvaluatedFunction2D[enrichmentFunctions.Count];
            for (int i = 0; i < enrichments.Length; ++i)
            {
                enrichments[i] = enrichmentFunctions[i].EvaluateAllAt(polarPoint, tipJacobians);
            }
            return enrichments;
        }


        // TODO: delete this
        public IReadOnlyList<CartesianPoint> IntersectionPointsForIntegration(XContinuumElement2D element)
        {
            throw new NotImplementedException();
        }
    }
}
