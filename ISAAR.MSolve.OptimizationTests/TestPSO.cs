using ISAAR.MSolve.Optimization;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.ParticleSwarmOptimization;
using ISAAR.MSolve.Optimization.Benchmarks.Mathematical;
using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Problems;
using System;

namespace ISAAR.MSolve.OptimizationTests
{
    class TestPSO
    {
        public static void Run()
        {
            OptimizationProblem optimizationProblem = new Ackley(2);

            ParticleSwarmOptimizationAlgorithm.Builder builder = new ParticleSwarmOptimizationAlgorithm.Builder(optimizationProblem);

            builder.SwarmSize = 10;
            builder.PhiP = 2.0;
            builder.PhiG = 2.0;
            builder.Omega = 0.2;
            builder.ConvergenceCriterion = new MaxFunctionEvaluations(10000);
            builder.Logger = new NoLogger();

            IOptimizationAlgorithm pso = builder.Build();
            IOptimizationAnalyzer analyzer = new OptimizationAnalyzer(pso);
            analyzer.Optimize();

            // Print results
            Console.WriteLine("\n Best Position:");
            for (int i = 0; i < optimizationProblem.Dimension; i++)
            {
                Console.WriteLine(String.Format(@"  x[{0}] = {1} ", i, pso.BestPosition[i]));
            }
            Console.WriteLine(String.Format(@"Best Fitness: {0}", pso.BestFitness));

            Console.Write("\nEnter any key to exit: ");
            Console.Read();
        }
    }
}
