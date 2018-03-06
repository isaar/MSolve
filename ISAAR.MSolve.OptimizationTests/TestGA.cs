using System;
using System.Linq;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling;
using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Problems;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations.Gaussian;
using ISAAR.MSolve.Optimization.Benchmarks.Mathematical;

namespace ISAAR.MSolve.OptimizationTests
{
    class TestGA
    {
        public static void Run()
        {
            //TestBinaryGA();
            TestRealCodedGA();
        }

        private static void TestBinaryGA()
        {
            // Define optim problem
            OptimizationProblem problem = new Rosenbrock();

            // Define optim algorithm
            var optimBuilder = new BinaryGeneticAlgorithmBuilder(problem);
            optimBuilder.Logger = new BestOfIterationLogger();

            // Setup logging
            //optimBuilder.Logger = new EmptyLogger();

            // Define convergence criteria 
            optimBuilder.ConvergenceCriterion = CompositeCriteria.OR(new MaxIterations(200), new MaxFunctionEvaluations(10000));

            // Define encoding
            optimBuilder.Encoding = new GrayCoding(problem, 16, 8);

            // Define population size
            optimBuilder.PopulationSize = 100;

            // Define selection strategy
            //optimBuilder.Selection = new RouletteWheelSelection<bool>(new InverseRankScaling<bool>(0.5));
            //optimBuilder.Selection = new TournamentSelection<bool>(2, false);
            //optimBuilder.Selection = new StochasticUniversalSampling<bool>(new InverseRankScaling<bool>(2.0
            //optimBuilder.Selection = new RemainderStochasticSamplingWithReplacement<bool>(new InverseRankScaling<bool>(1.0));
            optimBuilder.Selection = new TruncationSelection<bool>(0.5);

            // Define recombination strategy
            optimBuilder.Recombination = new SinglePointCrossover<bool>();

            // Define mutation strategy
            optimBuilder.Mutation = new BitFlipMutation(0.05);

            // Start optimization
            const int repetitions = 100;
            var solutions = new double[repetitions];
            for (int rep = 0; rep < repetitions; ++rep)
            {
                var optimAlgorithm = optimBuilder.BuildAlgorithm();
                optimAlgorithm.Solve();
                solutions[rep] = optimAlgorithm.BestFitness;
                Console.WriteLine("Best objective value: " + optimAlgorithm.BestFitness);
            }
            Console.WriteLine();
            Console.WriteLine("Average objective value: " + solutions.Average());
            Console.WriteLine();
            //Print results
            //Console.WriteLine("----------- History -----------");
            //optimBuilder.Logger.PrintToConsole();

            //Console.WriteLine("----------- Results -----------");
            //Console.WriteLine("Best objective value: " + optimAlgorithm.BestFitness);
            //Console.WriteLine("For continuous design variables: ");
            //PrintLineArray(optimAlgorithm.BestVariables.Item1);
            //Console.WriteLine("and integer design variables: ");
            //PrintLineArray(optimAlgorithm.BestVariables.Item2);
        }

        private static void TestRealCodedGA()
        {
            // Define optim problem
            OptimizationProblem problem = new Sphere();

            // Define optim algorithm
            var optimBuilder = new RealCodedGeneticAlgorithmBuilder(problem);
            optimBuilder.Logger = new BestOfIterationLogger();

            // Setup logging
            //optimBuilder.Logger = new EmptyLogger();

            // Define convergence criteria 
            optimBuilder.ConvergenceCriterion = CompositeCriteria.OR(new MaxIterations(200), new MaxFunctionEvaluations(100000));

            // Define population size
            optimBuilder.PopulationSize = 100;

            // Define selection strategy
            //optimBuilder.Selection = new RouletteWheelSelection<double>(new InverseRankScaling<double>(0.5));
            optimBuilder.Selection = new TournamentSelection<double>(2, false);
            //optimBuilder.Selection = new StochasticUniversalSampling<double>(new InverseRankScaling<double>(2.0
            //optimBuilder.Selection = new RemainderStochasticSamplingWithReplacement<double>(new InverseRankScaling<double>(1.0));
            //optimBuilder.Selection = new TruncationSelection<double>(0.5);

            // Define recombination strategy
            //optimBuilder.Recombination = new UniformCrossover<double>();
            optimBuilder.Recombination = new IntermediateCrossover();

            // Define mutation strategy
            //optimBuilder.Mutation = new ConstantGaussianMutation(problem, 1.0);
            optimBuilder.Mutation = new UniformMutation(problem, 0.01);

            // Start optimization
            const int repetitions = 100;
            var solutions = new double[repetitions];
            for (int rep = 0; rep < repetitions; ++rep)
            {
                var optimAlgorithm = optimBuilder.BuildAlgorithm();
                optimAlgorithm.Solve();
                solutions[rep] = optimAlgorithm.BestFitness;
                Console.WriteLine("Best objective value: " + optimAlgorithm.BestFitness);
            }
            Console.WriteLine();
            Console.WriteLine("Average objective value: " + solutions.Average());
            Console.WriteLine();
            //Print results
            //Console.WriteLine("----------- History -----------");
            //optimBuilder.Logger.PrintToConsole();

            //Console.WriteLine("----------- Results -----------");
            //Console.WriteLine("Best objective value: " + optimAlgorithm.BestFitness);
            //Console.WriteLine("For continuous design variables: ");
            //PrintLineArray(optimAlgorithm.BestVariables.Item1);
            //Console.WriteLine("and integer design variables: ");
            //PrintLineArray(optimAlgorithm.BestVariables.Item2);
        }
    }
}