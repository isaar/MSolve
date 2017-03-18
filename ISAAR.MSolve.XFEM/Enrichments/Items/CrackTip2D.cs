using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class CrackTip2D: AbstractEnrichmentItem2D
    {
        public ICartesianPoint2D TipCoordinates { get; private set; }
        public TipCoordinateSystem TipSystem { get; private set; }

        /// <summary>
        /// Angle (in rad) of local x to global x, in a counter clockwise rotation.
        /// </summary>
        public double LocalSystemOrientation { get; private set; }

        // The angle should be determined from the crack body curve, not by the user.
        public CrackTip2D(ICartesianPoint2D tipCoordinates, double localSystemOrientation)
        {
            this.TipCoordinates = tipCoordinates;
            this.LocalSystemOrientation = localSystemOrientation;
            this.TipSystem = new TipCoordinateSystem(this.TipCoordinates, this.LocalSystemOrientation);

            this.EnrichmentFunctions = new IEnrichmentFunction2D[]
            {
                new IsotropicBrittleTipFunctions2DAlternative.Func1(this),
                new IsotropicBrittleTipFunctions2DAlternative.Func2(this),
                new IsotropicBrittleTipFunctions2DAlternative.Func3(this),
                new IsotropicBrittleTipFunctions2DAlternative.Func4(this)
            };
            this.DOFs = new ArtificialDOFType[]
            {
                new ArtificialDOFType(EnrichmentFunctions[0], StandardDOFType.X),
                new ArtificialDOFType(EnrichmentFunctions[0], StandardDOFType.Y),
                new ArtificialDOFType(EnrichmentFunctions[1], StandardDOFType.X),
                new ArtificialDOFType(EnrichmentFunctions[1], StandardDOFType.Y),
                new ArtificialDOFType(EnrichmentFunctions[2], StandardDOFType.X),
                new ArtificialDOFType(EnrichmentFunctions[2], StandardDOFType.Y),
                new ArtificialDOFType(EnrichmentFunctions[3], StandardDOFType.X),
                new ArtificialDOFType(EnrichmentFunctions[3], StandardDOFType.Y),
            };
        }

        public override IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element)
        {
            throw new NotImplementedException("I should return the tip and the nodes of the element it is inside");
        }
    }
}
