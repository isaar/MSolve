using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;

// TODO: this class should not be associated with the whole crack geometry, just the part that stores the crack tip.
// Also it should be immutable, rather than having a settable TipSystem property. 
namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class CrackTipEnrichments2D : IEnrichmentItem2D
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
                new EnrichedDof(enrichmentFunctions[0], DisplacementDof.X),
                new EnrichedDof(enrichmentFunctions[0], DisplacementDof.Y),
                new EnrichedDof(enrichmentFunctions[1], DisplacementDof.X),
                new EnrichedDof(enrichmentFunctions[1], DisplacementDof.Y),
                new EnrichedDof(enrichmentFunctions[2], DisplacementDof.X),
                new EnrichedDof(enrichmentFunctions[2], DisplacementDof.Y),
                new EnrichedDof(enrichmentFunctions[3], DisplacementDof.X),
                new EnrichedDof(enrichmentFunctions[3], DisplacementDof.Y),
            };

        }

        public TipCoordinateSystem TipSystem { get; set; }

        public double[] EvaluateFunctionsAt(XNode2D node)
        {
            PolarPoint2D polarPoint = TipSystem.TransformPointGlobalCartesianToLocalPolar(node);
            var enrichments = new double[enrichmentFunctions.Count];
            for (int i = 0; i < enrichments.Length; ++i)
            {
                enrichments[i] = enrichmentFunctions[i].EvaluateAt(polarPoint);
            }
            return enrichments;
        }

        public EvaluatedFunction2D[] EvaluateAllAt(INaturalPoint2D point, XContinuumElement2D element,
             EvaluatedInterpolation2D interpolation)
        {
            PolarPoint2D polarPoint = TipSystem.TransformPointGlobalCartesianToLocalPolar(
                interpolation.TransformPointNaturalToGlobalCartesian(point));
            TipJacobians tipJacobians = TipSystem.CalculateJacobiansAt(polarPoint);

            var enrichments = new EvaluatedFunction2D[enrichmentFunctions.Count];
            for (int i = 0; i < enrichments.Length; ++i)
            {
                enrichments[i] = enrichmentFunctions[i].EvaluateAllAt(polarPoint, tipJacobians);
            }
            return enrichments;
        }


        // TODO: delete this
        public IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element)
        {
            throw new NotImplementedException();
        }
    }
}
