using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC;
using ISAAR.MSolve.Optimization.Logging;

//TODO: extend it to 3D.
//TODO: Should I store the number of elements and the prescribedVolume, instead of prescribedVolumeFraction?
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
            prescribedVolume = prescribedVolumeFraction * fem.CalVolume(Vector.CreateWithValue(fem.NumElements, 1.0));
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
            fem.AnalyzeModelWithDensities(densities.DoToAllEntries(xi => Math.Pow(xi, PenalizationExponent)));

            // Objective function and sensitivity analysis
            int numElements = fem.NumElements;
            int numLoadCases = fem.NumLoadCases;
            double compliance = 0.0;
            var sensitivities = Vector.CreateZero(numElements);
            for (int elem = 0; elem < numElements; ++elem)
            {
                IMatrixView stiffness = fem.GetBaseStiffnessOfElement(elem);
                double complianceCoeff = Math.Pow(densities[elem], PenalizationExponent);
                double sensitivityCoeff = -Math.Pow(densities[elem], PenalizationExponent - 1) * PenalizationExponent;
                for (int load = 0; load < numLoadCases; ++load)
                {
                    Vector displacements = fem.GetElementDisplacements(elem, load);
                    double unitCompliance = stiffness.Multiply(displacements).DotProduct(displacements);
                    compliance += complianceCoeff * unitCompliance;
                    sensitivities[elem] += sensitivityCoeff * unitCompliance;
                }
            }

            // Filter the densities and sensitivities
            //TODO: is it correct to modify the densities? They should be controlled by the optimization algorithm
            filter.FilterSensitivities(densities, ref sensitivities);

            if (Logger != null) Logger.Log(compliance);
            return (compliance, sensitivities);
        }

        private double EqualityConstraint(Vector densities) => fem.CalVolume(densities) - prescribedVolume; 
    }
}
