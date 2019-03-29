using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Tests;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Iterative;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using Xunit;

//TODO: add performance logging for solvers and gather all these in the same method.
namespace ISAAR.MSolve.Solvers.Tests
{
    public static class SingleSubdomainTests
    {
        [SkippableFact]
        internal static void TestDenseSolver()
        {
            Skip.IfNot(TestSettings.TestMkl, TestSettings.MessageWhenSkippingMKL);

            // Dense solver is too slow for a ~17.000 dof linear system, without MKL
            TestSettings.RunMultiproviderTest(LinearAlgebraProviderChoice.MKL, delegate ()
            {
                CantileverBeam benchmark = BuildCantileverBenchmark();

                var solverBuilder = new DenseMatrixSolver.Builder();
                DenseMatrixSolver solver = solverBuilder.BuildSolver(benchmark.Model);
                //solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()); // default

                RunAnalysisAndCheck(benchmark, solver);
            });
        }

        [Fact]
        internal static void TestPcgJacobiSolver()
        {
            CantileverBeam benchmark = BuildCantileverBenchmark();

            //LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
            var solverBuilder = new PcgSolver.Builder();
            //var pcgBuilder = new PcgAlgorithm.Builder();
            //solverBuilder.PcgAlgorithm = pcgBuilder.Build(); // default
            //solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()); // default
            PcgSolver solver = solverBuilder.BuildSolver(benchmark.Model);

            RunAnalysisAndCheck(benchmark, solver);
        }

        [Fact]
        internal static void TestPcgJacobiSolverWithAmdReordering()
        {
            CantileverBeam benchmark = BuildCantileverBenchmark();

            var solverBuilder = new PcgSolver.Builder();
            //var pcgBuilder = new PcgAlgorithm.Builder();
            //solverBuilder.PcgAlgorithm = pcgBuilder.Build();
            solverBuilder.DofOrderer = new DofOrderer(
                new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithCSparseAmd());
            PcgSolver solver = solverBuilder.BuildSolver(benchmark.Model);

            RunAnalysisAndCheck(benchmark, solver);
        }

        [Fact]
        internal static void TestSkylineSolver()
        {
            CantileverBeam benchmark = BuildCantileverBenchmark();

            var solverBuilder = new SkylineSolver.Builder();
            //solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()); // default
            ISolver_v2 solver = solverBuilder.BuildSolver(benchmark.Model);

            RunAnalysisAndCheck(benchmark, solver);
        }

        [Fact]
        internal static void TestSkylineSolverWithAmdReordering()
        {
            CantileverBeam benchmark = BuildCantileverBenchmark();

            var solverBuilder = new SkylineSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(
                new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithCSparseAmd());
            ISolver_v2 solver = solverBuilder.BuildSolver(benchmark.Model);
            RunAnalysisAndCheck(benchmark, solver);

        }

        [SkippableFact]
        internal static void TestSuiteSparseSolver()
        {
            Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

            CantileverBeam benchmark = BuildCantileverBenchmark();

            var solverBuilder = new SuiteSparseSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            using (SuiteSparseSolver solver = solverBuilder.BuildSolver(benchmark.Model))
            {
                RunAnalysisAndCheck(benchmark, solver);
            }
        }

        [SkippableFact]
        internal static void TestSuiteSparseSolverWithAmdReordering()
        {
            Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

            CantileverBeam benchmark = BuildCantileverBenchmark();

            var solverBuilder = new SuiteSparseSolver.Builder();
            //solverBuilder.DofOrderer = new DofOrderer(
            //    new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd()); // default
            using (SuiteSparseSolver solver = solverBuilder.BuildSolver(benchmark.Model))
            {
                RunAnalysisAndCheck(benchmark, solver);
            }
        }

        private static CantileverBeam BuildCantileverBenchmark()
        {
            var benchmarkBuilder = new CantileverBeam.Builder();
            //benchmarkBuilder.Length = 5.0;
            return benchmarkBuilder.BuildWithQuad4Elements(200, 10);
        }

        private static void RunAnalysisAndCheck(CantileverBeam benchmark, ISolver_v2 solver)
        {
            // Structural problem provider
            var provider = new ProblemStructural_v2(benchmark.Model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(benchmark.Model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(benchmark.Model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Check output
            double endDeflectionExpected = benchmark.CalculateEndDeflectionWithEulerBeamTheory();
            double endDeflectionComputed =
                benchmark.CalculateAverageEndDeflectionFromSolution(solver.LinearSystems.First().Value.Solution);
            var comparer = new ValueComparer(1E-2);
            Assert.True(comparer.AreEqual(endDeflectionExpected, endDeflectionComputed));
        }
    }
}
