using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.LinearAlgebra;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Materials.Crack;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Elements
{
    class XContinuumElementCrack2D: XContinuumElement2D
    {
        private readonly IIntegrationStrategy2D<XContinuumElementCrack2D> jIntegralStrategy;

        public XContinuumElementCrack2D(IsoparametricElementType2D type, IReadOnlyList<XNode2D> nodes,
            IIntegrationStrategy2D<XContinuumElement2D> integrationStrategy, 
            IIntegrationStrategy2D<XContinuumElementCrack2D> jIntegralStrategy): base(type, nodes, integrationStrategy)
        {
            this.jIntegralStrategy = jIntegralStrategy;
        }

        public void ComputeInteractionIntegrals(TipCoordinateSystem tipCoordinateSystem, CrackMaterial2D materialAtTip,
            double[] nodalDisplacementsX, double[] nodalDisplacementsY, double[] nodalWeights,
            out double integralMode1, out double integralMode2)
        {
            integralMode1 = 0.0;
            integralMode2 = 0.0;

            foreach (var entry in jIntegralStrategy.GetIntegrationPointsAndMaterials(this))
            {
                // Nomenclature: global = global cartesian system, natural = element natural system, 
                // local = tip local cartesian system  

                IFiniteElementMaterial2D material = entry.Value;
                GaussPoint2D naturalGP = entry.Key;
                EvaluatedInterpolation2D evaluatedInterpolation = Interpolation.EvaluateAt(Nodes, naturalGP);
                ICartesianPoint2D globalGP = evaluatedInterpolation.TransformPointNaturalToGlobalCartesian(naturalGP);

                // State 1
                DenseMatrix globalDisplacementGradState1 = CalculateDisplacementFieldGradient(evaluatedInterpolation, 
                    nodalDisplacementsX, nodalDisplacementsY);
                Tensor2D globalStressState1 = CalculateStressTensor(globalDisplacementGradState1, material);
                DenseMatrix localDisplacementGradState1 = tipCoordinateSystem.
                    TransformVectorFieldDerivativesGlobalCartesianToLocalCartesian(globalDisplacementGradState1);
                Tensor2D localStressTensorState1 = tipCoordinateSystem.
                    TransformTensorGlobalCartesianToLocalCartesian(globalStressState1);

                // Weight Function
                double[] globalWeightGradient = new double[2];
                for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
                {
                    var interpolationGradient = evaluatedInterpolation.GetGlobalCartesianDerivativesOf(Nodes[nodeIdx]);
                    globalWeightGradient[0] += interpolationGradient.Item1 * nodalWeights[nodeIdx];
                    globalWeightGradient[1] += interpolationGradient.Item2 * nodalWeights[nodeIdx];
                }
                double[] localWeightGradient = tipCoordinateSystem.
                    TransformScalarFieldDerivativesGlobalCartesianToLocalCartesian(globalWeightGradient);

                // State 2
                var auxiliaryStates = new AuxiliaryStates(materialAtTip, tipCoordinateSystem, globalGP);

                // Interaction integrals
                double integrandMode1 = ComputeJIntegrand(localWeightGradient, localDisplacementGradState1, 
                    localStressTensorState1, auxiliaryStates.DisplacementGradientMode1, 
                    auxiliaryStates.StrainTensorMode1, auxiliaryStates.StressTensorMode1);
                double integrandMode2 = ComputeJIntegrand(localWeightGradient, localDisplacementGradState1,
                    localStressTensorState1, auxiliaryStates.DisplacementGradientMode2,
                    auxiliaryStates.StrainTensorMode2, auxiliaryStates.StressTensorMode2);

                integralMode1 += integrandMode1 * evaluatedInterpolation.Jacobian.Determinant* naturalGP.Weight;
                integralMode2 += integrandMode2 * evaluatedInterpolation.Jacobian.Determinant * naturalGP.Weight;
            }
        }

        private static double ComputeJIntegrand(double[] weightGrad, DenseMatrix displGrad1, Tensor2D stress1,
            DenseMatrix displGrad2, Tensor2D strain2, Tensor2D stress2)
        {
            // Unrolled to greatly reduce mistakes. Alternatively Einstein notation products could be implementated
            // in Tensor2D, but still some parts would have to be unrolled.
            // Perhaps vector (and scalar) gradients should also be accessed by component and derivative variable.

            double strainEnergy = stress1.MultiplyColon(strain2);
            double parenthesis0 = stress1.XX * displGrad2[0, 0] + stress1.XY * displGrad2[1, 0]
                - stress2.XX * displGrad1[0, 0] - stress2.XY * displGrad1[1, 0] - strainEnergy;
            double parenthesis1 = stress1.XY * displGrad2[0, 0] + stress1.YY * displGrad2[1, 0]
                - stress2.XY * displGrad1[0, 0] - stress2.YY * displGrad1[1, 0];
            return parenthesis0 * weightGrad[0] + parenthesis1 * weightGrad[1];

        }
    }
}
