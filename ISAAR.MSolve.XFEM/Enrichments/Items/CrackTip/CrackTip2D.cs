using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Descriptions;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;


namespace ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip
{
    // TODO: there should be some linking between the crack tips and the crack body, outside sharing the same curve 
    // object.
    class CrackTip2D: AbstractEnrichmentItem2D
    {
        // TODO: The tip itself should not have to know where on the crack it is.
        // There must be a master crack class that knows where the tips, bodies and junctions are.
        public enum TipCurvePosition
        {
            CurveStart, CurveEnd, Both
        }

        private readonly IReadOnlyList<ITipFunction> enrichmentFunctions;
        private readonly TipCurvePosition tipPosition;
        private readonly IGeometryDescription2D discontinuity;
        private readonly ITipEnrichmentAreaStrategy enrichmentAreaStrategy;
        private readonly double magnificationOfJintegralRadius;
        private readonly IAuxiliaryStates auxiliaryStatesStrategy;
        private readonly ISIFCalculator sifCalculationStrategy;

        public IMesh2D<XNode2D, XContinuumElementCrack2D> Mesh { get; }

        // The next properties/fields need to be updated at each analysis step.
        public IReadOnlyList<XContinuumElementCrack2D> TipElements { get; private set; }
        public TipCoordinateSystem TipSystem { get; private set; }
        public ICartesianPoint2D TipCoordinates { get; private set; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tipPosition"></param>
        /// <param name="discontinuity"></param>
        /// <param name="enrichmentAreaStrategy"></param>
        /// <param name="magnificationOfJintegralRadius">The outer countour of the J-integral domain is defined as:
        ///     radius = magnification * sqrt(areaOfElementContainingTip). This parameter is the magnification. 
        ///     It should be at least 1.5 (see "Modeling quasi-static crack growth with the extended finite element 
        ///     method Part II: Numerical applications, Huang et al, 2003" page 7546). Usually values 2-3 are selected 
        ///     (see Ahmed thesis, 2009).</param>
        /// <param name="auxiliaryStatesStrategy"></param>
        /// <param name="sifCalculationStrategy"></param>
        public CrackTip2D(TipCurvePosition tipPosition, IGeometryDescription2D discontinuity,
            ITipEnrichmentAreaStrategy enrichmentAreaStrategy, double magnificationOfJintegralRadius,
            IAuxiliaryStates auxiliaryStatesStrategy, ISIFCalculator sifCalculationStrategy)
        {
            this.tipPosition = tipPosition;
            this.discontinuity = discontinuity;
            this.enrichmentAreaStrategy = enrichmentAreaStrategy;
            this.magnificationOfJintegralRadius = magnificationOfJintegralRadius; // TODO: Add checks for valid values
            this.auxiliaryStatesStrategy = auxiliaryStatesStrategy;
            this.sifCalculationStrategy = sifCalculationStrategy;

            this.enrichmentFunctions = new ITipFunction[]
            {
                new IsotropicBrittleTipFunctions2D.Func1(),
                new IsotropicBrittleTipFunctions2D.Func2(),
                new IsotropicBrittleTipFunctions2D.Func3(),
                new IsotropicBrittleTipFunctions2D.Func4()
            };
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

            UpdateTransform();
        }

        public override double[] EvaluateFunctionsAt(ICartesianPoint2D point)
        {
            PolarPoint2D polarPoint = TipSystem.TransformPointGlobalCartesianToLocalPolar(point);
            var enrichments = new double[enrichmentFunctions.Count];
            for (int i = 0; i < enrichments.Length; ++i)
            {
                enrichments[i] = enrichmentFunctions[i].EvaluateAt(polarPoint);
            }
            return enrichments;
        }

        public override EvaluatedFunction2D[] EvaluateAllAt(INaturalPoint2D point, IReadOnlyList<XNode2D> elementNodes,
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
            else if (tipPosition == TipCurvePosition.CurveStart) uniquePoints.Add(discontinuity.StartPoint);
            else throw new NotImplementedException("For now the tip can only be located at one end of the curve");

            return new List<ICartesianPoint2D>(uniquePoints);
        }

        private void UpdateTransform()
        {
            if (tipPosition == TipCurvePosition.CurveEnd)
            {
                TipCoordinates = discontinuity.EndPoint;
                TipSystem = new TipCoordinateSystem(TipCoordinates, discontinuity.EndPointOrientation());
            }
            else if (tipPosition == TipCurvePosition.CurveStart)
            {
                TipCoordinates = discontinuity.StartPoint;
                TipSystem = new TipCoordinateSystem(TipCoordinates, discontinuity.StartPointOrientation());
            }
            else
            {
                throw new NotImplementedException("For now the tip can only be located at one end of the curve");
            }
        }

        // TODO: what happens if the tip is on an element's edge/node?
        private void UpdateTipElements(ICartesianPoint2D newTip)
        {
            IReadOnlyList<XContinuumElementCrack2D> tipElements = Mesh.FindElementsContainingPoint(newTip, 
                TipElements[0]);
            if (tipElements.Count == 0) throw new NotImplementedException("New tip is outside of domain");
            else TipElements = tipElements;
        }

        private void ComputeSIFS(Model2D model, double[] totalFreeDisplacements, double[] totalConstrainedDisplacements,
            out double sifMode1, out double sifMode2)
        {
            double interactionIntegralMode1 = 0.0, interactionIntegralMode2 = 0.0;
            foreach(var pair in FindJintegralElementsAndNodalWeights())
            {
                XContinuumElementCrack2D element = pair.Key;
                double[] nodalWeights = pair.Value;
                double[] elementDisplacements = model.DofEnumerator.ExtractDisplacementVectorOfElementFromGlobal(
                    element, totalFreeDisplacements, totalConstrainedDisplacements);

                double partialIntegralMode1, partialIntegralMode2;
                element.ComputeInteractionIntegrals(TipSystem, auxiliaryStatesStrategy, elementDisplacements,
                    nodalWeights, out partialIntegralMode1, out partialIntegralMode2);

                interactionIntegralMode1 += partialIntegralMode1;
                interactionIntegralMode2 += partialIntegralMode2;
            }

            sifMode1 = sifCalculationStrategy.CalculateSIF(interactionIntegralMode1);
            sifMode2 = sifCalculationStrategy.CalculateSIF(interactionIntegralMode2);
        }

        private IReadOnlyDictionary<XContinuumElementCrack2D, double[]> FindJintegralElementsAndNodalWeights()
        {
            Circle2D outerContour = 
                new Circle2D(TipCoordinates, ComputeRadiusOfJintegralOuterContour());
            IReadOnlyList<XContinuumElementCrack2D> intersectedElements = 
                Mesh.FindElementsIntersectedByCircle(outerContour, TipElements[0]);

            var elementsAndWeights = new Dictionary<XContinuumElementCrack2D, double[]>();
            foreach (var element in intersectedElements)
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
                    else // Node lies inside or exactly on the circle
                    {
                        nodalWeights[nodeIdx] = 1.0;
                    }
                }
                elementsAndWeights.Add(element, nodalWeights);
            }
            return elementsAndWeights;
        }

        private double ComputeRadiusOfJintegralOuterContour()
        {
            double maxTipElementArea = -1.0;
            foreach (var element in TipElements)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes);
                double elementArea = outline.ComputeArea();
                if (elementArea > maxTipElementArea) maxTipElementArea = elementArea;
            }
            return magnificationOfJintegralRadius * Math.Sqrt(maxTipElementArea);
        }
    }
}
