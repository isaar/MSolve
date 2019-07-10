using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Analysis;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Filtering;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.MaterialInterpolation;

//TODO: This probably works for 3D as it is. If not, extend it.
//TODO: Should I store the prescribedVolume, instead of prescribedVolumeFraction?
//TODO: Add a builder.
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP
{
    public class TopologySimpLinear2D
    {
        private readonly ILinearFemAnalysis fem;
        private readonly IDensityFilter filter;
        private readonly IMaterialInterpolation materialInterpolation;
        private readonly OptimalityCriteriaBuilder optimAlgorithmBuilder;
        private readonly double prescribedVolumeFraction;
        private OptimalityCriteria optimAlgorithm; //TODO: extend this
        private double prescribedVolume;

        public TopologySimpLinear2D(ILinearFemAnalysis fem, OptimalityCriteriaBuilder optimAlgorithmBuilder,
            IDensityFilter filter, IMaterialInterpolation materialInterpolation, double prescribedVolumeFraction)
        {
            this.fem = fem;
            this.optimAlgorithmBuilder = optimAlgorithmBuilder;
            this.filter = filter;
            this.materialInterpolation = materialInterpolation;
            this.prescribedVolumeFraction = prescribedVolumeFraction;
        }

        public ObjectiveFunctionLogger Logger { get; set; } //TODO: Extend this

        public void Initialize()
        {
            fem.Initialize();
            prescribedVolume = prescribedVolumeFraction * fem.CalculateTotalVolume(Vector.CreateWithValue(fem.NumElements, 1.0));
            optimAlgorithm = optimAlgorithmBuilder.BuildOptimizer(ObjectiveFunction, EqualityConstraint, 
                materialInterpolation.MinDensity, materialInterpolation.MaxDensity);
        }

        public (double minCompliance, Vector bestDensities) Optimize()
        {
            (Vector bestDensities, double minCompliance) = 
                optimAlgorithm.Optimize(Vector.CreateWithValue(fem.NumElements, prescribedVolumeFraction));
            return (minCompliance, bestDensities);
        }

        private (double compliance, Vector sensitivities) ObjectiveFunction(Vector densities)
        {
            // FEM analysis. Use the material property (e.g Young modulus) for the current element densities.
            IVectorView materialProperties = densities.DoToAllEntries(materialInterpolation.CalcMaterialProperty);
            fem.UpdateModelAndAnalyze(materialProperties);

            // Objective function and sensitivity analysis
            int numElements = fem.NumElements;
            int numLoadCases = fem.NumLoadCases;
            double compliance = 0.0;
            var sensitivities = Vector.CreateZero(numElements);
            for (int e = 0; e < numElements; ++e)
            {
                IMatrixView stiffness = fem.GetElementStiffnessForUnitMaterial(e);
                double complianceCoeff = materialProperties[e];
                double sensitivityCoeff = - materialInterpolation.CalcMaterialPropertyDerivative(densities[e]);
                for (int load = 0; load < numLoadCases; ++load)
                {
                    // If K0(e) is the element's stiffness matrix for E = 1 and E(x(e)) is the Young modulus for its current 
                    // density:
                    // Compliance: c += E(x(e)) * ( U(e)^T * K0(e) * U(e) )
                    // Sensitivity: dc(e) -= dE(x(e))/dx(e) * (  U(e)^T * K0(e) * U(e) ),
                    Vector displacements = fem.GetElementDisplacements(e, load);
                    double unitCompliance = stiffness.Multiply(displacements).DotProduct(displacements);
                    compliance += complianceCoeff * unitCompliance;
                    sensitivities[e] += sensitivityCoeff * unitCompliance;
                }
            }

            // Filter the densities and sensitivities
            //TODO: is it correct to modify the densities? They should be controlled by the optimization algorithm
            filter.FilterSensitivities(densities, ref sensitivities);

            if (Logger != null) Logger.Log(compliance);
            return (compliance, sensitivities);
        }

        private double EqualityConstraint(Vector densities) => fem.CalculateTotalVolume(densities) - prescribedVolume;
    }
}
