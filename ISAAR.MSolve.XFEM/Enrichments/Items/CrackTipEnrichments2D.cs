using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.Descriptions;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class CrackTipEnrichments2D : IEnrichmentItem2D
    {
        private readonly IReadOnlyList<ITipFunction> enrichmentFunctions;
        public IReadOnlyList<ArtificialDOFType> DOFs { get; }
        public ICrackDescription crackDescription;

        public CrackTipEnrichments2D(ICrackDescription crackDescription): this(crackDescription, new ITipFunction[] 
            {
                new IsotropicBrittleTipFunctions2D.Func1(),
                new IsotropicBrittleTipFunctions2D.Func2(),
                new IsotropicBrittleTipFunctions2D.Func3(),
                new IsotropicBrittleTipFunctions2D.Func4()
            })
        {
        }

        public CrackTipEnrichments2D(ICrackDescription crackDescription, IReadOnlyList<ITipFunction> enrichmentFunctions)
        {
            this.crackDescription = crackDescription;
            this.enrichmentFunctions = enrichmentFunctions;
            this.DOFs = new ArtificialDOFType[]
            {
                new ArtificialDOFType(enrichmentFunctions[0], StandardDOFType.X),
                new ArtificialDOFType(enrichmentFunctions[0], StandardDOFType.Y),
                new ArtificialDOFType(enrichmentFunctions[1], StandardDOFType.X),
                new ArtificialDOFType(enrichmentFunctions[1], StandardDOFType.Y),
                new ArtificialDOFType(enrichmentFunctions[2], StandardDOFType.X),
                new ArtificialDOFType(enrichmentFunctions[2], StandardDOFType.Y),
                new ArtificialDOFType(enrichmentFunctions[3], StandardDOFType.X),
                new ArtificialDOFType(enrichmentFunctions[3], StandardDOFType.Y),
            };

        }

        public double[] EvaluateFunctionsAt(XNode2D node)
        {
            PolarPoint2D polarPoint = crackDescription.TipSystem.TransformPointGlobalCartesianToLocalPolar(node);
            var enrichments = new double[enrichmentFunctions.Count];
            for (int i = 0; i < enrichments.Length; ++i)
            {
                enrichments[i] = enrichmentFunctions[i].EvaluateAt(polarPoint);
            }
            return enrichments;
        }

        public EvaluatedFunction2D[] EvaluateAllAt(INaturalPoint2D point, IReadOnlyList<XNode2D> elementNodes,
             EvaluatedInterpolation2D interpolation)
        {
            PolarPoint2D polarPoint = crackDescription.TipSystem.TransformPointGlobalCartesianToLocalPolar(
                interpolation.TransformPointNaturalToGlobalCartesian(point));
            TipJacobians tipJacobians = crackDescription.TipSystem.CalculateJacobiansAt(polarPoint);

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
