using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Solvers.MenkBordas;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class SlopeSolvers
    {
        private delegate ISolver CreateSolver(Slope benchmark);

        public static void Run()
        {
            TestSolver();
            //CompareSolvers();
        }

        private static void TestSolver()
        {
            Slope.Builder builder = SetupBenchmark();
            builder.LsmOutputDirectory = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Slope\Results";
            Slope benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();

            // Solvers
            var solver = CreateCholeskySuiteSparseSolver(benchmark);
            //var solver = CreateCholeskyAMDSolver(benchmark);
            //var solver = CreateMenkBordasSolver(benchmark);
            IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solver);

            //Reanalysis solvers
            //IReadOnlyList<ICartesianPoint2D> crackPath;
            //using (var solver = CreateReanalysisRebuildingSolver(benchmark))
            //using (var solver = CreateReanalysisSolver(benchmark))
            //{
            //    crackPath = benchmark.Analyze(solver);
            //}

            Console.WriteLine("Crack path:");
            foreach (var point in crackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
        }

        private static void CompareSolvers()
        {
            int numRepetitions = 10;

            // Define the benchmark problem
            Slope.Builder builder = SetupBenchmark();

            // Define the solvers
            string[] solvers =
            {
                "SuiteSparse Cholesky", "Reanalysis Pais"
            };
            CreateSolver[] callbacks =
            {
                CreateCholeskySuiteSparseSolver, CreateReanalysisRebuildingSolver
            };
            var totalTimes = new long[solvers.Length];

            // Run the analyses by alternating them to reduce bias.
            for (int t = 0; t < numRepetitions; ++t)
            {
                Console.WriteLine($"Repetition: {t}");
                for (int s = 0; s < solvers.Length; ++s)
                {
                    Slope benchmark = builder.BuildBenchmark();
                    benchmark.InitializeModel();
                    var solver = callbacks[s](benchmark);
                    benchmark.Analyze(solver);
                    long totalTime = solver.Logger.CalcTotalTime();
                    Console.WriteLine($"Solver {solvers[s]}: total time = {totalTime} ms.");
                    totalTimes[s] += totalTime;
                }
                Console.WriteLine();
            }

            // Print statistics
            for (int s = 0; s < solvers.Length; ++s)
            {
                Console.WriteLine($"Solver {solvers[s]}: Total time = {totalTimes[s] / numRepetitions} ms");
            }
        }

        public static Slope.Builder SetupBenchmark()
        {
            string meshPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Slope\Meshes\slope.msh";
            double growthLength = 2.0; // Must be sufficientrly larger than the element size.
            var builder = new Slope.Builder(meshPath, growthLength);
            builder.MaxIterations = 10;
            

            // TODO: enter the fixed propagator here, perhaps by solving the benchmark once.
            // These are for mesh: ...
            //var propagationLogger = new PropagationLogger();
            //propagationLogger.GrowthAngles.Add(-0.0349434289780521);
            //propagationLogger.GrowthLengths.Add(0.3);
            //propagationLogger.GrowthAngles.Add(-0.0729848767244629);
            //propagationLogger.GrowthLengths.Add(0.3);
            //propagationLogger.GrowthAngles.Add(-0.125892740180586);
            //propagationLogger.GrowthLengths.Add(0.3);
            //propagationLogger.GrowthAngles.Add(-0.200116860828933);
            //propagationLogger.GrowthLengths.Add(0.3);
            //propagationLogger.GrowthAngles.Add(-0.258299391791769);
            //propagationLogger.GrowthLengths.Add(0.3);
            //propagationLogger.GrowthAngles.Add(-0.264803465603906);
            //propagationLogger.GrowthLengths.Add(0.3);
            //propagationLogger.GrowthAngles.Add(-0.201411670680886);
            //propagationLogger.GrowthLengths.Add(0.3);
            //propagationLogger.GrowthAngles.Add(-0.123234163953279);
            //propagationLogger.GrowthLengths.Add(0.3);
            //propagationLogger.GrowthAngles.Add(-0.0816346096256186);
            //propagationLogger.GrowthLengths.Add(0.3);
            //builder.KnownPropagation = propagationLogger;

            return builder;
        }

        private static ISolver CreateMenkBordasSolver(Slope benchmark)
        {
            return new MenkBordasSolver(benchmark.Model, benchmark.Decomposer, 1000000, double.Epsilon);
        }

        private static ISolver CreateCholeskySuiteSparseSolver(Slope benchmark)
        {
            return new CholeskySuiteSparseSolver(benchmark.Model);
        }

        private static ISolver CreateCholeskyAMDSolver(Slope benchmark)
        {
            return new CholeskyAMDSolver(benchmark.Model);
        }

        private static ReanalysisRebuildingSolver CreateReanalysisRebuildingSolver(Slope benchmark)
        {
            return new ReanalysisRebuildingSolver(benchmark.Model, benchmark.Crack);
        }

        private static ReanalysisSolver CreateReanalysisSolver(Slope benchmark)
        {
            return new ReanalysisSolver(benchmark.Model, benchmark.Crack);
        }
    }
}
