using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class CrackBodyEnrichment2D : IEnrichmentItem2D
    {
        private readonly IHeavisideFunction2D enrichmentFunction;
        public IReadOnlyList<ArtificialDOFType> DOFs { get; }
        public IExteriorCrack crackDescription;

        public CrackBodyEnrichment2D(IExteriorCrack crackDescription): this(crackDescription, new SignFunction2D())
        {
        }

        public CrackBodyEnrichment2D(IExteriorCrack crackDescription, IHeavisideFunction2D enrichmentFunction)
        {
            this.crackDescription = crackDescription;
            this.enrichmentFunction = enrichmentFunction;
            this.DOFs = new ArtificialDOFType[] {
                new ArtificialDOFType(enrichmentFunction, StandardDOFType.X),
                new ArtificialDOFType(enrichmentFunction, StandardDOFType.Y)
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
