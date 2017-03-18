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
using ISAAR.MSolve.XFEM.Geometry.Descriptions;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    // TODO: there should be some linking between the crack tips and the crack body, outside sharing the same curve 
    // object.
    class CrackTip2D: AbstractEnrichmentItem2D
    {
        
        private readonly TipPosition tipPosition;
        private readonly IGeometryDescription2D discontinuity;

        // This needs to be updated at each analysis step.
        public TipCoordinateSystem TipSystem { get; private set; }

        /// <summary>
        /// Angle (in rad) of local x to global x, in a counter clockwise rotation.
        /// </summary>
        public double LocalSystemOrientation { get; private set; }

        // The angle should be determined from the crack body curve, not by the user.
        public CrackTip2D(TipPosition tipPosition, IGeometryDescription2D discontinuity)
        {
            this.tipPosition = tipPosition;
            this.discontinuity = discontinuity;

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

            UpdateTransform();
        }

        public void UpdateTransform()
        {
            if (tipPosition == TipPosition.CurveEnd)
            {
                this.TipSystem = new TipCoordinateSystem(discontinuity.EndPoint, discontinuity.EndPointOrientation());
            }
            else
            {
                throw new NotImplementedException("For now the tip can only be located at the curve's end");
            }
        }

        // TODO: this needs reworking. E.g. making sure there are no identical points might already be 
        // done by the geometry class
        public override IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element)
        {
            var uniquePoints = new HashSet<ICartesianPoint2D>(element.Nodes);
            
            foreach (ICartesianPoint2D point in discontinuity.IntersectionWith(element))
            {
                uniquePoints.Add(point);
            }

            if (tipPosition == TipPosition.CurveEnd) uniquePoints.Add(discontinuity.EndPoint);
            else throw new NotImplementedException("For now the tip can only be located at the curve's end");

            return new List<ICartesianPoint2D>(uniquePoints);
        }

        // TODO: a more polymorhic design would be better
        public enum TipPosition
        {
            CurveStart, CurveEnd
        }
    }
}
