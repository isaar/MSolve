using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Solvers.MenkBordas;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class DCBSolvers
    {
        private const string outputPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_DCB\Plots";

        private delegate ISolver CreateSolver(DCB benchmark);

        public static void Run()
        {
            TestSolver();
            //CompareSolvers();
        }

        private static void TestSolver()
        {
            DCB.Builder builder = SetupBenchmark();
            builder.LsmOutputDirectory = outputPath;
            DCB benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();

            // Solvers
            //var solver = CreateSkylineSolver(benchmark);
            //var solver = CreateCholeskySuiteSparseSolver(benchmark);
            //var solver = CreateNoReanalysisSolver(benchmark);
            //var solver = CreateCholeskyAMDSolver(benchmark);
            //var solver = CreatePCGSolver(benchmark);
            //var solver = CreateMinresSolver(benchmark);
            var solver = CreateMenkBordasSolver(benchmark);
            IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solver);

            //Reanalysis solvers
            //IReadOnlyList<ICartesianPoint2D> crackPath;
            ////using (var solver = CreateReanalysisRebuildingSolver(benchmark))
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
            DCB.Builder builder = SetupBenchmark();

            // Define the solvers
            string[] solvers =
                {
                "Skyline", "Jacobi Preconditioned CG", "SuiteSparse Cholesky", "Reanalysis Pais"
            };
            CreateSolver[] callbacks =
            {
                CreateSkylineSolver, CreatePCGSolver, CreateCholeskySuiteSparseSolver, CreateReanalysisRebuildingSolver
            };
            //var totalTimes = new long[solvers.Length];

            // Run the analyses by alternating them to reduce bias.
            for (int t = 0; t < numRepetitions; ++t)
            {
                Console.WriteLine($"Repetition: {t}");
                for (int s = 0; s < solvers.Length; ++s)
                {
                    DCB benchmark = builder.BuildBenchmark();
                    benchmark.InitializeModel();
                    var solver = callbacks[s](benchmark);
                    benchmark.Analyze(solver);
                    //long totalTime = solver.Logger.CalcTotalTime();
                    //Console.WriteLine($"Solver {solvers[s]}: total time = {totalTime} ms.");
                    //totalTimes[s] += totalTime;
                }
                Console.WriteLine();
            }

            // Print statistics
            //for (int s = 0; s < solvers.Length; ++s)
            //{
            //    Console.WriteLine($"Solver {solvers[s]}: Total time = {totalTimes[s] / numRepetitions} ms");
            //}
        }

        public static DCB.Builder SetupBenchmark()
        {
            double growthLength = 0.3; //0.3, 0.4 is good
            double elementSize = 0.08; //0.08, 0.2 is good
            IMeshProvider meshProvider = new DCBUniformMeshProvider(elementSize);
            //double fineElementSize = 0.11;
            //IMeshProvider meshProvider = new DCBRefinedMeshProvider(fineElementSize, 10 * fineElementSize);
            //IMeshProvider meshProvider = new GmshMeshProvider(@"C: \Users\Serafeim\Desktop\GMSH\dcb.msh");
            var builder = new DCB.Builder(growthLength, meshProvider);
            builder.UseLSM = true;
            //TODO: fix a bug that happens when the crack has almost reached the boundary, is inside but no tip elements are 
            //      found. It happens at iteration 10.
            builder.MaxIterations = 10;
            //builder.KnownPropagation = null; // TODO: enter the fixed propagator here, perhaps by solving the benchmark once.

            // These are for uniform mesh, element size = 0.08;
            var propagationLogger = new PropagationLogger();
            propagationLogger.GrowthAngles.Add(-0.0349434289780521);
            propagationLogger.GrowthLengths.Add(0.3);
            propagationLogger.GrowthAngles.Add(-0.0729848767244629);
            propagationLogger.GrowthLengths.Add(0.3);
            propagationLogger.GrowthAngles.Add(-0.125892740180586);
            propagationLogger.GrowthLengths.Add(0.3);
            propagationLogger.GrowthAngles.Add(-0.200116860828933);
            propagationLogger.GrowthLengths.Add(0.3);
            propagationLogger.GrowthAngles.Add(-0.258299391791769);
            propagationLogger.GrowthLengths.Add(0.3);
            propagationLogger.GrowthAngles.Add(-0.264803465603906);
            propagationLogger.GrowthLengths.Add(0.3);
            propagationLogger.GrowthAngles.Add(-0.201411670680886);
            propagationLogger.GrowthLengths.Add(0.3);
            propagationLogger.GrowthAngles.Add(-0.123234163953279);
            propagationLogger.GrowthLengths.Add(0.3);
            propagationLogger.GrowthAngles.Add(-0.0816346096256186);
            propagationLogger.GrowthLengths.Add(0.3);
            //builder.KnownPropagation = propagationLogger;

            return builder;
        }

        private static ISolver CreateSkylineSolver(DCB benchmark)
        {
            return new SkylineSolver(benchmark.Model);
        }

        private static ISolver CreatePCGSolver(DCB benchmark)
        {
            return new PCGSolver(benchmark.Model, 1, 1e-10);
        }

        private static ISolver CreateMinresSolver(DCB benchmark)
        {
            return new MinresSolver(benchmark.Model, 1, 1e-10, 0);
        }

        private static ISolver CreateMenkBordasSolver(DCB benchmark)
        {
            return new MenkBordasSolver(benchmark.Model, benchmark.Crack, benchmark.Decomposer, 1000000, 1e-10, 
                new StandardPreconditionerCholesky.Builder(benchmark.Model), new EnrichedPreconditioningNaive(), outputPath);
        }

        private static ISolver CreateCholeskySuiteSparseSolver(DCB benchmark)
        {
            return new CholeskySuiteSparseSolver(benchmark.Model);
        }

        private static ISolver CreateCholeskyAMDSolver(DCB benchmark)
        {
            return new CholeskyAMDSolver(benchmark.Model);
        }

        private static ReanalysisRebuildingSolver CreateReanalysisRebuildingSolver(DCB benchmark)
        {
            return new ReanalysisRebuildingSolver(benchmark.Model, benchmark.Crack);
        }

        private static ReanalysisSolver CreateReanalysisSolver(DCB benchmark)
        {
            return new ReanalysisSolver(benchmark.Model, benchmark.Crack);
        }

        private static NoReanalysisSolver CreateNoReanalysisSolver(DCB benchmark)
        {
            return new NoReanalysisSolver(benchmark.Model, benchmark.Crack);
        }


    }
}
