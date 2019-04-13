using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.DifferentialEvolution;
using ISAAR.MSolve.Optimization.Benchmarks.Unconstrained;
using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Problems;
using Xunit;

namespace ISAAR.MSolve.Optimization.Tests
{
    public static class TestDE
    {
        [Fact]
        public static void Run()
        {
            int seed = 1;
            var rng = new Random(seed);

            OptimizationProblem optimizationProblem = new Rosenbrock();

            var builder = new DifferentialEvolutionAlgorithm.Builder(optimizationProblem);
            builder.PopulationSize = 100;
            builder.MutationFactor = 0.6;
            builder.CrossoverProbability = 0.9;
            builder.ConvergenceCriterion = new MaxFunctionEvaluations(100000);
            builder.RandomNumberGenerator = rng;
            IOptimizationAlgorithm de = builder.Build();

            IOptimizationAnalyzer analyzer = new OptimizationAnalyzer(de);
            analyzer.Optimize();

            double expectedFitness = 0.0;
            var expectedDesign = Vector.CreateWithValue(optimizationProblem.Dimension, 1.0);
            Assert.Equal(expectedFitness, de.BestFitness, 10);
            Assert.True(Vector.CreateFromArray(de.BestPosition).Equals(expectedDesign, 1E-6));
        }
    }
}
