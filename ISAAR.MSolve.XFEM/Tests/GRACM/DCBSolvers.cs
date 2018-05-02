using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Solvers;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class DCBSolvers
    {
        private delegate ISolver CreateSolver(DCB benchmark);

        public static void Main()
        {
            TestSolver();
            //CompareSolvers();
        }

        private static void TestSolver()
        {
            DCB.Builder builder = SetupBenchmark();
            builder.UseLSM = true; // Explicit crack results in a singular matrix. I probably broke sth. TODO: solve this
            DCB benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();

            // Skyline
            //var solver = CreateSkylineSolver(benchmark);
            //IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solver);

            // SuiteSparse Cholesky
            //var solver = CreateCholeskySuiteSparseSolver(benchmark);
            //IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solver);

            // Reanalysis
            //IReadOnlyList<ICartesianPoint2D> crackPath;
            //using (var solver = CreateReanalysisSolver(benchmark))
            //{
            //    crackPath = benchmark.Analyze(solver);
            //}

            // Benchmark AMD
            var solver = CreateNoReanalysisSolver(benchmark);
            IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solver);


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
                CreateSkylineSolver, CreatePCGSolver, CreateCholeskySuiteSparseSolver, CreateReanalysisSolver
            };
            var totalTimes = new long[solvers.Length];

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

        private static DCB.Builder SetupBenchmark()
        {
            double growthLength = 0.3;
            double fineElementSize = 0.08;
            IMeshProvider meshProvider = new DCBRefinedMeshProvider(fineElementSize, 10 * fineElementSize);
            //IMeshProvider meshProvider = new GmshMeshProvider(@"C: \Users\Serafeim\Desktop\GMSH\dcb.msh");
            var builder = new DCB.Builder(growthLength, meshProvider);
            builder.UseLSM = true;
            //TODO: fix a bug that happens when the crack has almost reached the boundary, is inside but no tip elements are 
            //      found. It happens at iteration 10.
            builder.MaxIterations = 10;
            builder.KnownPropagation = null; // TODO: enter the fixed propagator here, perhaps by solving the benchmark once.

            return builder;
        }

        private static ISolver CreateSkylineSolver(DCB benchmark)
        {
            return new SkylineSolver(benchmark.Model);
        }

        private static ISolver CreatePCGSolver(DCB benchmark)
        {
            return new PCGSolver(benchmark.Model, 1, 1e-8);
        }

        private static ISolver CreateCholeskySuiteSparseSolver(DCB benchmark)
        {
            return new CholeskySuiteSparseSolver(benchmark.Model);
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
