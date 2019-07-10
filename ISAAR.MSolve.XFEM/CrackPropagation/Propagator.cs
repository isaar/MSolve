using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.CrackPropagation
{
    public class Propagator: IPropagator
    {
        private readonly IMesh2D<XNode, XContinuumElement2D> mesh;
        private readonly double magnificationOfJintegralRadius;
        private readonly IAuxiliaryStates auxiliaryStatesStrategy;
        private readonly ISIFCalculator sifCalculationStrategy;
        private readonly ICrackGrowthDirectionLaw2D growthDirectionLaw;
        private readonly ICrackGrowthLengthLaw2D growthLengthLaw;

        public PropagationLogger Logger { get; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="crack"></param>
        /// <param name="magnificationOfJintegralRadius">The outer countour of the J-integral domain is defined as:
        ///     radius = magnification * sqrt(areaOfElementContainingTip). This parameter is the magnification. 
        ///     It should be at least 1.5 (see "Modeling quasi-static crack growth with the extended finite element 
        ///     method Part II: Numerical applications, Huang et al, 2003" page 7546). Usually values 2-3 are selected 
        ///     (see Ahmed thesis, 2009).</param>
        /// <param name="auxiliaryStatesStrategy"></param>
        /// <param name="sifCalculationStrategy"></param>
        public Propagator(IMesh2D<XNode, XContinuumElement2D> mesh, double magnificationOfJintegralRadius,
            IAuxiliaryStates auxiliaryStatesStrategy, ISIFCalculator sifCalculationStrategy,
            ICrackGrowthDirectionLaw2D growthDirectionLaw, ICrackGrowthLengthLaw2D growthLengthLaw)
        {
            this.mesh = mesh;
            this.magnificationOfJintegralRadius = magnificationOfJintegralRadius;
            this.auxiliaryStatesStrategy = auxiliaryStatesStrategy;
            this.sifCalculationStrategy = sifCalculationStrategy;
            this.growthDirectionLaw = growthDirectionLaw;
            this.growthLengthLaw = growthLengthLaw;
            this.Logger = new PropagationLogger();
        }

        public (double growthAngle, double growthLength) Propagate(Dictionary<int, Vector> totalFreeDisplacements, 
            CartesianPoint crackTip, TipCoordinateSystem tipSystem, IReadOnlyList<XContinuumElement2D> tipElements)
        {
            // TODO: Also check if the sifs do not violate the material toughness
            (double sifMode1, double sifMode2) = ComputeSIFS(totalFreeDisplacements, crackTip, tipSystem, tipElements);
            double growthAngle = growthDirectionLaw.ComputeGrowthAngle(sifMode1, sifMode2);
            double growthLength = growthLengthLaw.ComputeGrowthLength(sifMode1, sifMode2);
            Logger.GrowthAngles.Add(growthAngle);
            Logger.GrowthLengths.Add(growthLength);
            return (growthAngle, growthLength);

        }

        private (double sifMode1, double sifMode2) ComputeSIFS(Dictionary<int, Vector> totalFreeDisplacements, 
            CartesianPoint crackTip, TipCoordinateSystem tipSystem, IReadOnlyList<XContinuumElement2D> tipElements)
        {
            double interactionIntegralMode1 = 0.0, interactionIntegralMode2 = 0.0;
            IReadOnlyDictionary<XContinuumElement2D, double[]> elementWeights = 
                FindJintegralElementsAndNodalWeights(crackTip, tipSystem, tipElements);
            foreach (var pair in elementWeights)
            {
                XContinuumElement2D element = pair.Key;
                double[] nodalWeights = pair.Value;

                //TODO: This needs refactoring ASAP.
                XSubdomain subdomain = element.Subdomain;
                double[] elementDisplacements = 
                    subdomain.CalculateElementDisplacements(element, totalFreeDisplacements[subdomain.ID]);
                (double[] standardElementDisplacements, double[] enrichedElementDisplacements) = 
                    element.SeparateStdEnrVector(elementDisplacements);

                //Vector standardElementDisplacements = dofOrderer.ExtractDisplacementVectorOfElementFromGlobal(
                //    element, totalFreeDisplacements, totalConstrainedDisplacements);
                //Vector enrichedElementDisplacements = dofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(
                //    element, totalFreeDisplacements);

                double partialIntegralMode1, partialIntegralMode2;
                ComputeInteractionIntegrals(element, Vector.CreateFromArray(standardElementDisplacements), 
                    Vector.CreateFromArray(enrichedElementDisplacements), nodalWeights, tipSystem, 
                    out partialIntegralMode1, out partialIntegralMode2);

                interactionIntegralMode1 += partialIntegralMode1;
                interactionIntegralMode2 += partialIntegralMode2;
            }

            double sifMode1 = sifCalculationStrategy.CalculateSIF(interactionIntegralMode1);
            double sifMode2 = sifCalculationStrategy.CalculateSIF(interactionIntegralMode2);

            Logger.InteractionIntegralsMode1.Add(interactionIntegralMode1);
            Logger.InteractionIntegralsMode2.Add(interactionIntegralMode2);
            Logger.SIFsMode1.Add(sifMode1);
            Logger.SIFsMode2.Add(sifMode2);

            return (sifMode1, sifMode2);
        }

        private IReadOnlyDictionary<XContinuumElement2D, double[]> FindJintegralElementsAndNodalWeights(
            CartesianPoint crackTip, TipCoordinateSystem tipSystem, IReadOnlyList<XContinuumElement2D> tipElements)
        {
            Circle2D outerContour = 
                new Circle2D(crackTip, ComputeRadiusOfJintegralOuterContour(tipSystem, tipElements));
            IReadOnlyList<XContinuumElement2D> intersectedElements =
                mesh.FindElementsIntersectedByCircle(outerContour, tipElements[0]);

            var elementsAndWeights = new Dictionary<XContinuumElement2D, double[]>();
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
                    CirclePointPosition pos = outerContour.FindRelativePositionOfPoint((CartesianPoint)element.Nodes[nodeIdx]);
                    if (pos == CirclePointPosition.Outside)
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

        // TODO: This method should directly return the elements and take care of cases near the domain boundaries (see Ahmed)
        // TODO: The J-integral radius should not exceed the last crack segment's length
        public double ComputeRadiusOfJintegralOuterContour(TipCoordinateSystem tipSystem, 
            IReadOnlyList<XContinuumElement2D> tipElements)
        {
            double maxTipElementArea = -1.0;
            foreach (var element in tipElements)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes.Select(node => (CartesianPoint)node).ToArray());
                double elementArea = outline.ComputeArea();
                if (elementArea > maxTipElementArea) maxTipElementArea = elementArea;
            }
            return magnificationOfJintegralRadius * Math.Sqrt(maxTipElementArea);
        }

        private void ComputeInteractionIntegrals(XContinuumElement2D element, Vector standardNodalDisplacements,
            Vector enrichedNodalDisplacements, double[] nodalWeights, TipCoordinateSystem tipSystem,
            out double integralMode1, out double integralMode2)
        {
            integralMode1 = 0.0;
            integralMode2 = 0.0;
            foreach (GaussPoint naturalGP in element.JintegralStrategy.GenerateIntegrationPoints(element))
            {
                // Nomenclature: global = global cartesian system, natural = element natural system, 
                // local = tip local cartesian system  
                
                EvalInterpolation2D evaluatedInterpolation = 
                    element.Interpolation.EvaluateAllAt(element.Nodes, naturalGP);
                CartesianPoint globalGP = evaluatedInterpolation.TransformPointNaturalToGlobalCartesian();
                Matrix constitutive = 
                    element.Material.CalculateConstitutiveMatrixAt(naturalGP, evaluatedInterpolation);

                // State 1
                Matrix2by2 globalDisplacementGradState1 = element.CalculateDisplacementFieldGradient(
                    naturalGP, evaluatedInterpolation, standardNodalDisplacements, enrichedNodalDisplacements);
                Tensor2D globalStressState1 = element.CalculateStressTensor(globalDisplacementGradState1, constitutive);
                Matrix2by2 localDisplacementGradState1 = tipSystem.
                    TransformVectorFieldDerivativesGlobalCartesianToLocalCartesian(globalDisplacementGradState1);
                Tensor2D localStressTensorState1 = tipSystem.
                    TransformTensorGlobalCartesianToLocalCartesian(globalStressState1);

                // Weight Function
                // TODO: There should be a method InterpolateScalarGradient(double[] nodalValues) in EvaluatedInterpolation
                // TODO: Rewrite this as a shapeGradients (matrix) * nodalWeights (vector) operation.
                var globalWeightGradient = Vector2.CreateZero();
                for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
                {
                    globalWeightGradient.AxpyIntoThis(
                        evaluatedInterpolation.ShapeGradientsCartesian.GetRow(nodeIdx), // Previously: GetGlobalCartesianDerivativesOf(element.Nodes[nodeIdx])
                        nodalWeights[nodeIdx]);
                }
                Vector2 localWeightGradient = tipSystem.
                    TransformScalarFieldDerivativesGlobalCartesianToLocalCartesian(globalWeightGradient);

                // State 2
                // TODO: XContinuumElement shouldn't have to pass tipCoordinate system to auxiliaryStates. 
                // It would be better to have CrackTip handle this and the coordinate transformations. That would also 
                // obey LoD, but a lot of wrapper methods would be required.
                AuxiliaryStatesTensors auxiliary = auxiliaryStatesStrategy.ComputeTensorsAt(globalGP, tipSystem);

                // Interaction integrals
                double integrandMode1 = ComputeJIntegrand(localWeightGradient, localDisplacementGradState1,
                    localStressTensorState1, auxiliary.DisplacementGradientMode1,
                    auxiliary.StrainTensorMode1, auxiliary.StressTensorMode1);
                double integrandMode2 = ComputeJIntegrand(localWeightGradient, localDisplacementGradState1,
                    localStressTensorState1, auxiliary.DisplacementGradientMode2,
                    auxiliary.StrainTensorMode2, auxiliary.StressTensorMode2);

                integralMode1 += integrandMode1 * evaluatedInterpolation.Jacobian.DirectDeterminant * naturalGP.Weight;
                integralMode2 += integrandMode2 * evaluatedInterpolation.Jacobian.DirectDeterminant * naturalGP.Weight;
            }
        }

        private static double ComputeJIntegrand(Vector2 weightGrad, Matrix2by2 displGrad1, Tensor2D stress1,
            Matrix2by2 displGrad2, Tensor2D strain2, Tensor2D stress2)
        {
            // Unrolled to greatly reduce mistakes. Alternatively Einstein notation products could be implementated
            // in Tensor2D (like the tensor-tensor multiplication is), but still some parts would have to be unrolled.
            // Perhaps vector (and scalar) gradients should also be accessed by component and derivative variable.

            double strainEnergy = stress1.MultiplyColon(strain2);
            double parenthesis0 = stress1.XX * displGrad2[0, 0] + stress1.XY * displGrad2[1, 0]
                + stress2.XX * displGrad1[0, 0] + stress2.XY * displGrad1[1, 0] - strainEnergy;
            double parenthesis1 = stress1.XY * displGrad2[0, 0] + stress1.YY * displGrad2[1, 0]
                + stress2.XY * displGrad1[0, 0] + stress2.YY * displGrad1[1, 0];
            return parenthesis0 * weightGrad[0] + parenthesis1 * weightGrad[1];
        }
    }
}
