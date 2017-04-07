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
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Descriptions;
//using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;


namespace ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip
{
    // TODO: there should be some linking between the crack tips and the crack body, outside sharing the same curve 
    // object.
    class CrackTip2D: AbstractEnrichmentItem2D
    {
        // TODO: a more polymorhic design would be better
        public enum TipCurvePosition
        {
            CurveStart, CurveEnd
        }

        private readonly TipCurvePosition tipPosition;
        private readonly IGeometryDescription2D discontinuity;

        // The next properties/fields need to be updated at each analysis step.
        public TipCoordinateSystem TipSystem { get; private set; }
        public ICartesianPoint2D TipCoordinates { get; private set; }

        /// <summary>
        /// Angle (in rad) of local x to global x, in a counter clockwise rotation.
        /// </summary>
        public double LocalSystemOrientation { get; private set; }

        private XContinuumElementCrack2D TipElement { get; set; }

        // The angle should be determined from the crack body curve, not by the user.
        public CrackTip2D(TipCurvePosition tipPosition, IGeometryDescription2D discontinuity)
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
            if (tipPosition == TipCurvePosition.CurveEnd)
            {
                TipCoordinates = discontinuity.EndPoint;
                TipSystem = new TipCoordinateSystem(TipCoordinates, discontinuity.EndPointOrientation());
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

            if (tipPosition == TipCurvePosition.CurveEnd) uniquePoints.Add(discontinuity.EndPoint);
            else throw new NotImplementedException("For now the tip can only be located at the curve's end");

            return new List<ICartesianPoint2D>(uniquePoints);
        }

        private double ComputeJintegralOuterRadius()
        {
            throw new NotImplementedException("Needs input from the constructor and is different for fixed enrichment");
        }

        private IReadOnlyDictionary<XContinuumElementCrack2D, double[]> FindJintegralElementsAndNodalWeights(
            Model2D model)
        {
            Circle2D outerContour = new Circle2D(TipCoordinates, ComputeJintegralOuterRadius());

            var elementsAndWeights = new Dictionary<XContinuumElementCrack2D, double[]>();
            foreach (Element2D element in model.Elements) // TODO: Reduce the O(elementsCount) complexity by using a better mesh class.
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes);
                if (outline.IntersectsWithCircle(outerContour))
                {
                    // The relative position of the circle and the nodes was already calculated when checking the
                    // circle-element intersection, but that method should be decoupled from assigning the nodal 
                    // weights, even at the cost of some duplicate operations. What could be done more efficiently is 
                    // caching the nodes and weights already processed by previous elements, but even then the cost of
                    // processing each node will be increased by the lookup.

                    double[] nodalWeights = new double[element.Nodes.Count];
                    for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
                    {
                        if (outerContour.FindRelativePositionOfPoint(element.Nodes[nodeIdx]) == CirclePointPosition.Outside)
                        {
                            nodalWeights[nodeIdx] = 0.0;
                        }
                        else // Node is inside or exactly on the circle
                        {
                            nodalWeights[nodeIdx] = 1.0;
                        }
                    }
                    elementsAndWeights.Add((XContinuumElementCrack2D)(element.ElementType), nodalWeights); // TODO: This is horrible. Model should be generic.
                }
            }
            return elementsAndWeights;
        }
    }
}
