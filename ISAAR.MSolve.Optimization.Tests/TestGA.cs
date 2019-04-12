using System;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;
using ISAAR.MSolve.Optimization.Benchmarks.Unconstrained;
using ISAAR.MSolve.Optimization.Commons;
using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Initialization;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Problems;
using Troschuetz.Random.Generators;
using Xunit;

namespace ISAAR.MSolve.Optimization.Tests
{
    public static class TestGA
    {
        //[Fact] // This algorithm is too slow for Rosenbrock test. It does converge though.
        public static void TestBinaryGA()
        {
            int seed = 1;
            var rng = new StandardGenerator(seed);
            //var rng = RandomNumberGenerationUtilities.troschuetzRandom;

            // Define optim problem
            OptimizationProblem problem = new Rosenbrock();

            // Define optim algorithm
            var optimBuilder = new BinaryGeneticAlgorithmBuilder(problem);
            optimBuilder.Logger = new BestOfIterationLogger();

            // Setup logging
            //optimBuilder.Logger = new EmptyLogger();

            // Define convergence criteria 
            optimBuilder.ConvergenceCriterion = new MaxFunctionEvaluations(1000000);

            // Define encoding
            optimBuilder.Encoding = new GrayCoding(problem, 16, 8);

            // Initialization
            optimBuilder.Initializer = new RealUniformRandomInitializer(problem, rng);

            // Define population size
            optimBuilder.PopulationSize = 100;

            // Define selection strategy
            //optimBuilder.Selection = new RouletteWheelSelection<bool>(new InverseRankScaling<bool>(0.5));
            //optimBuilder.Selection = new TournamentSelection<bool>(2, false);
            //optimBuilder.Selection = new StochasticUniversalSampling<bool>(new InverseRankScaling<bool>(2.0
            //optimBuilder.Selection = new RemainderStochasticSamplingWithReplacement<bool>(new InverseRankScaling<bool>(1.0));
            optimBuilder.Selection = new TruncationSelection<bool>(rng, 0.5);

            // Define recombination strategy
            optimBuilder.Recombination = new SinglePointCrossover<bool>(rng);

            // Define mutation strategy
            optimBuilder.Mutation = new BitFlipMutation(0.05, rng);

            // Start optimization
            GeneticAlgorithm<bool> optimAlgorithm = optimBuilder.BuildAlgorithm();
            optimAlgorithm.Solve();

            // Check results
            double expectedFitness = 0.0;
            var expectedDesign = Vector.CreateWithValue(problem.Dimension, 1.0);
            Assert.Equal(expectedFitness, optimAlgorithm.BestFitness, 10);
            Assert.True(Vector.CreateFromArray(optimAlgorithm.BestPosition).Equals(expectedDesign, 1E-6));
        }

        [Fact]
        public static void TestRealCodedGA()
        {
            int seed = 2;
            var rng = new StandardGenerator(seed);
            RandomNumberGenerationUtilities.troschuetzRandom = rng; // hack to make sure the same seed is used everwhere in the GA. TODO: Find out which part of the GA messes up.

            // Define optim problem
            OptimizationProblem problem = new Sphere();

            // Define optim algorithm
            var optimBuilder = new RealCodedGeneticAlgorithmBuilder(problem);
            optimBuilder.Logger = new BestOfIterationLogger();

            // Setup logging
            //optimBuilder.Logger = new EmptyLogger();

            // Initialization
            optimBuilder.Initializer = new RealUniformRandomInitializer(problem, rng);

            // Define convergence criteria 
            optimBuilder.ConvergenceCriterion = CompositeCriteria.OR(new MaxIterations(200), new MaxFunctionEvaluations(1000000));

            // Define population size
            optimBuilder.PopulationSize = 100;

            // Define selection strategy
            //optimBuilder.Selection = new RouletteWheelSelection<double>(new InverseRankScaling<double>(0.5));
            optimBuilder.Selection = new TournamentSelection<double>(rng, 2, false);
            //optimBuilder.Selection = new StochasticUniversalSampling<double>(new InverseRankScaling<double>(2.0
            //optimBuilder.Selection = new RemainderStochasticSamplingWithReplacement<double>(new InverseRankScaling<double>(1.0));
            //optimBuilder.Selection = new TruncationSelection<double>(0.5);

            // Define recombination strategy
            //optimBuilder.Recombination = new UniformCrossover<double>();
            optimBuilder.Recombination = new IntermediateCrossover(rng);

            // Define mutation strategy
            //optimBuilder.Mutation = new ConstantGaussianMutation(problem, 1.0);
            optimBuilder.Mutation = new UniformMutation(problem, 0.01, rng);

            // Start optimization
            GeneticAlgorithm<double> optimAlgorithm = optimBuilder.BuildAlgorithm();
            optimAlgorithm.Solve();

            // Check results
            double expectedFitness = 0.0;
            var expectedDesign = Vector.CreateWithValue(problem.Dimension, 0.0);
            Assert.Equal(expectedFitness, optimAlgorithm.BestFitness, 2);
            Assert.True(Vector.CreateFromArray(optimAlgorithm.BestPosition).Equals(expectedDesign, 1E-1));
        }
    }
}