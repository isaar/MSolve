using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Analysis;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Filters;

//TODO: This probably works for 3D as it is. If not, extend it.
//TODO: Should I store the prescribedVolume, instead of prescribedVolumeFraction?
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP
{
    public class TopologySimpLinear2D
    {
        private const double maxDensity = 1.0;

        private readonly ILinearFemAnalysis fem;
        private readonly IDensityFilter filter;
        private readonly double prescribedVolumeFraction;
        private OptimalityCriteria optimAlgorithm; //TODO: extend this
        private double prescribedVolume;

        public TopologySimpLinear2D(ILinearFemAnalysis fem, IDensityFilter filter, double prescribedVolumeFraction)
        {
            this.fem = fem;
            this.filter = filter;
            this.prescribedVolumeFraction = prescribedVolumeFraction;
        }

        public ObjectiveFunctionLogger Logger { get; set; } //TODO: Extend this

        //TODO: perhaps a builder for these default values
        public double MinDensity { get; set; } = 0.001;
        public double PenalizationExponent { get; set; } = 3.0;

        public void Initialize()
        {
            fem.Initialize();
            prescribedVolume = prescribedVolumeFraction * fem.CalculateTotalVolume(Vector.CreateWithValue(fem.NumElements, 1.0));
            optimAlgorithm = new OptimalityCriteria(ObjectiveFunction, EqualityConstraint, MinDensity, maxDensity);
        }

        public (double minCompliance, Vector bestDensities) Optimize()
        {
            (double minCompliance, Vector bestDensities) = 
                optimAlgorithm.Optimize(Vector.CreateWithValue(fem.NumElements, prescribedVolumeFraction));
            return (minCompliance, bestDensities);
        }

        private (double compliance, Vector sensitivities) ObjectiveFunction(Vector densities)
        {
            // FEM analysis. Use the penalized densities.
            Vector penalizedDensities = densities.DoToAllEntries(xi => Math.Pow(xi, PenalizationExponent));
            fem.AnalyzeModelWithDensities(penalizedDensities);

            // Objective function and sensitivity analysis
            int numElements = fem.NumElements;
            int numLoadCases = fem.NumLoadCases;
            double compliance = 0.0;
            var sensitivities = Vector.CreateZero(numElements);
            for (int e = 0; e < numElements; ++e)
            {
                IMatrixView stiffness = fem.GetPenalizedStiffnessOfElement(e, penalizedDensities);
                double sensitivityCoeff = -PenalizationExponent * Math.Pow(densities[e], PenalizationExponent - 1) 
                    / penalizedDensities[e];
                for (int load = 0; load < numLoadCases; ++load)
                {
                    // In theory, if x(e) is the density of an element, p is the penalty and K(e) is its stiffness matrix before 
                    // penalizing its material properties:
                    // Compliance: c += x(e)^p * ( U(e)^T * K(e) * U(e) )
                    // Sensitivity: dc(e) -= p * x(e)^(p-1) * (  U(e)^T * K(e) * U(e) ),
                    // However, calculating both K(e) here and x(e)^p*K(e) during global matrix assembly is cumbersome.
                    // Therefore we will only calculate x(e)^p*K(e) and rewrite the equations above as:
                    // Compliance: c += U(e)^T * x(e)^p*K(e) * U(e)
                    // Sensitivity: dc(e) -= p * ( x(e)^(p-1) / x(e)^p ) * (  U(e)^T * x(e)^p*K(e) * U(e) ),
                    Vector displacements = fem.GetElementDisplacements(e, load);
                    double elementCompliance = stiffness.Multiply(displacements).DotProduct(displacements);
                    compliance += elementCompliance;
                    sensitivities[e] += sensitivityCoeff * elementCompliance;
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
