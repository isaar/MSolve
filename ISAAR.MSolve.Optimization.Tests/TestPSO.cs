using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.ParticleSwarmOptimization;
using ISAAR.MSolve.Optimization.Benchmarks.Unconstrained;
using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Problems;
using Xunit;

namespace ISAAR.MSolve.Optimization.Tests
{
    public static class TestPSO
    {
        [Fact]
        public static void Run()
        {
            int seed = 1;
            var rng = new Random(seed);

            OptimizationProblem optimizationProblem = new Ackley(2);

            var builder = new ParticleSwarmOptimizationAlgorithm.Builder(optimizationProblem);
            builder.SwarmSize = 10;
            builder.PhiP = 2.0;
            builder.PhiG = 2.0;
            builder.Omega = 0.2;
            builder.ConvergenceCriterion = new MaxFunctionEvaluations(10000);
            builder.Logger = new NoLogger();
            builder.RandomNumberGenerator = rng;

            IOptimizationAlgorithm pso = builder.Build();
            IOptimizationAnalyzer analyzer = new OptimizationAnalyzer(pso);
            analyzer.Optimize();

            double expectedFitness = 0.0;
            var expectedDesign = Vector.CreateZero(optimizationProblem.Dimension);
            Assert.Equal(expectedFitness, pso.BestFitness, 3);
            Assert.True(Vector.CreateFromArray(pso.BestPosition).Equals(expectedDesign, 1E-4));
        }
    }
}
