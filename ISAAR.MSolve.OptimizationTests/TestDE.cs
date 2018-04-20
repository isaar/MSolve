using System;
using ISAAR.MSolve.Optimization;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.DifferentialEvolution;
using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Benchmarks.Mathematical;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.OptimizationTests
{
    public class TestDE
    {
        public static void Run()
        {
            OptimizationProblem optimizationProblem = new Rosenbrock();

            DifferentialEvolutionAlgorithm.Builder builder = new DifferentialEvolutionAlgorithm.Builder(optimizationProblem);
            builder.PopulationSize = 100;
            builder.MutationFactor = 0.6;
            builder.CrossoverProbability = 0.9;
            builder.ConvergenceCriterion = new MaxFunctionEvaluations(100000);
            IOptimizationAlgorithm de = builder.Build();

            IOptimizationAnalyzer analyzer = new OptimizationAnalyzer(de);
            analyzer.Optimize();

            // Print results
            Console.WriteLine("\n Best Position:");
            for (int i = 0; i < optimizationProblem.Dimension; i++)
            {
                Console.WriteLine(String.Format(@"  x[{0}] = {1} ", i, de.BestPosition[i]));
            }
            Console.WriteLine(String.Format(@"Best Fitness: {0}", de.BestFitness));
        }
    }
}
