using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.DifferentialEvolution;
using ISAAR.MSolve.Optimization.Benchmarks.Mathematical.Constrained;
using ISAAR.MSolve.Optimization.Constraints.Penalties;
using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.OptimizationTests
{
    public class TestDEConstrained
    {
        public static void Run()
        {
            OptimizationProblem optimizationProblem = new S_CRES();

            var builder = new DifferentialEvolutionAlgorithmConstrained.Builder(optimizationProblem);
            builder.PopulationSize = 20;
            builder.MutationFactor = 0.6;
            builder.CrossoverProbability = 0.9;
            builder.ConvergenceCriterion = new MaxFunctionEvaluations(100000);
            builder.Penalty = new DeathPenalty();
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
