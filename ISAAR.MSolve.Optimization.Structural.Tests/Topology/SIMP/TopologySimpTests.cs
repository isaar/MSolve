using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP;
using Xunit;

namespace ISAAR.MSolve.Optimization.Structural.Tests.Topology.SIMP
{
    /// <summary>
    /// Tests for <see cref="TopologySimpLinear2D"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class TopologySimpTests
    {
        [Fact]
        private static void TestCantileverBeam()
        {
            // Parameters
            int numElementsX = 32, numElementsY = 20;
            double filterAreaRadius = 1.2, volumeFraction = 0.4, penalty = 3.0;

            // Define the optimization
            var material = new ElasticMaterial2D(StressState2D.PlaneStress) { YoungModulus = 1.0, PoissonRatio = 0.3 };
            var fem = new LinearFemAnalysisUniform2D(numElementsX, numElementsY, material,
                LinearFemAnalysisUniform2D.BoundaryConditions.ShortCantilever);
            var filter = new MeshIndependentSensitivityFilter2DUniform(numElementsX, numElementsY, filterAreaRadius);
            var simp = new TopologySimpLinear2D(fem, filter, volumeFraction);
            simp.PenalizationExponent = penalty;
            var logger = new ObjectiveFunctionLogger();
            simp.Logger = logger;

            // Run the optimization
            simp.Initialize();
            (double compliance, Vector densities) = simp.Optimize();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(TopologySimp99LinesTests.ComplianceHistoryForCantilever);
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-11));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_shortcantilever_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Vector densitiesExpected = reader.ReadFile(densitiesPath).Reshape(false);
            Assert.True(densitiesExpected.Equals(densities, 1E-10));
        }

        [Fact]
        private static void TestCantileverBeamWith2LoadCases()
        {
            // Parameters
            int numElementsX = 30, numElementsY = 30;
            double filterAreaRadius = 1.2, volumeFraction = 0.4, penalty = 3.0;

            // Define the optimization
            var material = new ElasticMaterial2D(StressState2D.PlaneStress) { YoungModulus = 1.0, PoissonRatio = 0.3 };
            var fem = new LinearFemAnalysisUniform2D(numElementsX, numElementsY, material,
                LinearFemAnalysisUniform2D.BoundaryConditions.Cantilever2LoadCases);
            var filter = new MeshIndependentSensitivityFilter2DUniform(numElementsX, numElementsY, filterAreaRadius);
            var simp = new TopologySimpLinear2D(fem, filter, volumeFraction);
            simp.PenalizationExponent = penalty;
            var logger = new ObjectiveFunctionLogger();
            simp.Logger = logger;

            // Run the optimization
            simp.Initialize();
            (double compliance, Vector densities) = simp.Optimize();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(
                TopologySimp99LinesTests.ComplianceHistoryForCantileverWith2LoadCases);
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-11));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_cantilever_loadcases_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Vector densitiesExpected = reader.ReadFile(densitiesPath).Reshape(false); 
            Assert.True(densitiesExpected.Equals(densities, 1E-11));
        }

        [Fact]
        private static void TestMbbBeam()
        {
            // Parameters
            int numElementsX = 60, numElementsY = 20;
            double filterAreaRadius = 1.5, volumeFraction = 0.5, penalty = 3.0;

            // Define the optimization
            var material = new ElasticMaterial2D(StressState2D.PlaneStress) { YoungModulus = 1.0, PoissonRatio = 0.3 };
            var fem = new LinearFemAnalysisUniform2D(numElementsX, numElementsY, material, 
                LinearFemAnalysisUniform2D.BoundaryConditions.MbbBeam);
            var filter = new MeshIndependentSensitivityFilter2DUniform(numElementsX, numElementsY, filterAreaRadius);
            var simp = new TopologySimpLinear2D(fem, filter, volumeFraction);
            simp.PenalizationExponent = penalty;
            var logger = new ObjectiveFunctionLogger();
            simp.Logger = logger;

            // Run the optimization
            simp.Initialize();
            (double compliance, Vector densities) = simp.Optimize();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(TopologySimp99LinesTests.ComplianceHistoryForMbbBeam);
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-10));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_MBBbeam_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Vector densitiesExpected = reader.ReadFile(densitiesPath).Reshape(false);
            Assert.True(densitiesExpected.Equals(densities, 1E-9));
        }
    }
}
