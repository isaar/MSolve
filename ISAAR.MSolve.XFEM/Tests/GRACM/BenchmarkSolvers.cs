using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Solvers.MenkBordas;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class BenchmarkSolvers
    {
        private delegate ISolver CreateSolver(IBenchmark benchmark);

        public static void Run()
        {
            SingleTest();
            //BenchmarkSolver();
        }


        private static void SingleTest()
        {
            IBenchmarkBuilder builder = Fillet.SetupBenchmark();
            IBenchmark benchmark = builder.BuildBenchmark();
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

        private static void BenchmarkSolver()
        {

        }

        private static ISolver CreateMenkBordasSolver(IBenchmark benchmark)
        {
            return new MenkBordasSolver(benchmark.Model, benchmark.Decomposer, 1000000, double.Epsilon);
        }

        private static ISolver CreateCholeskySuiteSparseSolver(IBenchmark benchmark)
        {
            return new CholeskySuiteSparseSolver(benchmark.Model);
        }

        private static ISolver CreateCholeskyAMDSolver(IBenchmark benchmark)
        {
            return new CholeskyAMDSolver(benchmark.Model);
        }

        private static ReanalysisRebuildingSolver CreateReanalysisRebuildingSolver(IBenchmark benchmark)
        {
            return new ReanalysisRebuildingSolver(benchmark.Model, benchmark.EnrichedArea, benchmark.Crack);
        }

        private static ReanalysisSolver CreateReanalysisSolver(IBenchmark benchmark)
        {
            return new ReanalysisSolver(benchmark.Model, benchmark.EnrichedArea, benchmark.Crack);
        }
    }
}
