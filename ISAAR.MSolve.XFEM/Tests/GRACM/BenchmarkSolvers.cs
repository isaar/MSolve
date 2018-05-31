using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.CrackGeometry;
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
            //IBenchmarkBuilder builder = Fillet.SetupBenchmark(true, true);
            IBenchmarkBuilder builder = Holes.SetupBenchmark(false, true);


            IBenchmark benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();

            // Solvers
            //var solver = CreateCholeskySuiteSparseSolver(benchmark);
            var solver = CreateMenkBordasSolver(benchmark);
            //var solver = CreatePCGSolver(benchmark);
            //var solver = CreateCholeskyAMDSolver(benchmark);
            benchmark.Analyze(solver);

            //Reanalysis solvers
            //var solver = CreateReanalysisRebuildingSolver(benchmark);
            //var solver = CreateReanalysisSolver(benchmark);
            //benchmark.Analyze(solver);
            //solver.Dispose();

            // Timing output path
            string timingOutputPath = builder.TimingOutputDirectory + "\\" + solver.Logger.SolverName + "_results.txt";
            solver.Logger.WriteSumsToFile(timingOutputPath, benchmark.Name, true);
        }

        private static void BenchmarkSolver()
        {
            int numRepetitions = 10;

            /// Define the benchmark problem once
            //IBenchmarkBuilder builder = Fillet.SetupBenchmark(false, false);
            IBenchmarkBuilder builder = Holes.SetupBenchmark(false, false);

            /// Choose solver
            //CreateSolver solverFunc = CreateCholeskyAMDSolver;
            //CreateSolver solverFunc = CreateReanalysisSolver;
            CreateSolver solverFunc = CreatePCGSolver;
            //CreateSolver solverFunc = CreateMenkBordasSolver;

            /// Call once to load all necessary DLLs
            IBenchmark firstTry = builder.BuildBenchmark();
            firstTry.InitializeModel();
            ISolver firstSolver = solverFunc(firstTry);
            firstTry.Analyze(firstSolver);

            /// Timing output path
            string timingOutputPath = builder.TimingOutputDirectory + "\\" + firstSolver.Logger.SolverName + "_results.txt";

            /// Actually time the solver
            for (int t = 0; t < numRepetitions; ++t)
            {
                /// Create a new benchmark at each iteration
                IBenchmark benchmark = builder.BuildBenchmark();
                benchmark.InitializeModel();

                /// Run the analysis
                ISolver solver = solverFunc(benchmark);
                benchmark.Analyze(solver);

                /// Write the timing results
                solver.Logger.WriteSumsToFile(timingOutputPath, benchmark.Name, true);

                /// Dispose any unmanaged memory
                if (solver is IDisposable handle) handle.Dispose();
            }
        }

        private static ISolver CreateCholeskySuiteSparseSolver(IBenchmark benchmark)
        {
            return new CholeskySuiteSparseSolver(benchmark.Model);
        }

        private static ISolver CreateCholeskyAMDSolver(IBenchmark benchmark)
        {
            return new CholeskyAMDSolver(benchmark.Model);
        }

        private static ISolver CreateMenkBordasSolver(IBenchmark benchmark)
        {
            return new MenkBordasSolver(benchmark.Model, benchmark.Crack, benchmark.Decomposer, 1000000, 1e-10,
                benchmark.PlotDirectory);
        }

        private static ISolver CreatePCGSolver(IBenchmark benchmark)
        {
            return new PCGSolver(benchmark.Model, 1, 1e-10);
        }

        private static ReanalysisRebuildingSolver CreateReanalysisRebuildingSolver(IBenchmark benchmark)
        {
            return new ReanalysisRebuildingSolver(benchmark.Model, benchmark.Crack, benchmark.PossibleEnrichments);
        }

        private static ReanalysisSolver CreateReanalysisSolver(IBenchmark benchmark)
        {
            return new ReanalysisSolver(benchmark.Model, benchmark.Crack, benchmark.PossibleEnrichments);
        }
    }
}
