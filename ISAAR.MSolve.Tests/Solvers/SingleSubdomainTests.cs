﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Dense;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.PCG;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;

//TODO: add performance logging for solvers and gather all these in the same method.
namespace ISAAR.MSolve.Tests.Solvers
{
    public static class SingleSubdomainTests
    {
        [Fact]
        internal static void TestDenseSolver()
        {
            var solverBuilder = new DenseMatrixSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            RunCantileverBenchmark(model => solverBuilder.BuildSolver(model));
        }

        [Fact]
        internal static void TestPcgJacobiSolver()
        {
            var pcgBuilder = new PcgAlgorithm.Builder();
            var solverBuilder = new PcgSolver.Builder(pcgBuilder.Build());
            solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            RunCantileverBenchmark(model => solverBuilder.BuildSolver(model));
        }

        [Fact]
        internal static void TestPcgJacobiSolverWithAmdReordering()
        {
            var pcgBuilder = new PcgAlgorithm.Builder();
            var solverBuilder = new PcgSolver.Builder(pcgBuilder.Build());
            solverBuilder.DofOrderer = new DofOrderer(
                new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
            RunCantileverBenchmark(model => solverBuilder.BuildSolver(model));
        }

        [Fact]
        internal static void TestSkylineSolver()
        {
            var solverBuilder = new SkylineSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            RunCantileverBenchmark(model => solverBuilder.BuildSolver(model));
        }

        [Fact]
        internal static void TestSkylineSolverWithAmdReordering()
        {
            var solverBuilder = new SkylineSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(
                new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
            RunCantileverBenchmark(model => solverBuilder.BuildSolver(model));
        }

        [Fact]
        internal static void TestSuiteSparseSolver()
        {
            var solverBuilder = new SuiteSparseSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            RunCantileverBenchmark(model => solverBuilder.BuildSolver(model));
        }

        [Fact]
        internal static void TestSuiteSparseSolverWithAmdReordering()
        {
            var solverBuilder = new SuiteSparseSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(
                new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
            RunCantileverBenchmark(model => solverBuilder.BuildSolver(model));
        }

        private static void RunCantileverBenchmark(Func<Model_v2, ISolver_v2> buildSolver)
        {
            var benchmarkBuilder = new CantileverBeam.Builder();
            benchmarkBuilder.Length = 5.0;
            CantileverBeam benchmark = benchmarkBuilder.BuildWithQuad4Elements(160, 8);

            // Solver
            ISolver_v2 solver = buildSolver(benchmark.Model);

            // Structural problem provider
            var provider = new ProblemStructural_v2(benchmark.Model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(solver);
            var parentAnalyzer = new StaticAnalyzer_v2(benchmark.Model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Check output
            double endDeflectionExpected = benchmark.CalculateEndDeflectionWithEulerBeamTheory();
            double endDeflectionComputed = 
                benchmark.CalculateAverageEndDeflectionFromSolution(solver.LinearSystems.First().Solution);
            var comparer = new ValueComparer(1E-2);
            //Assert.True(comparer.AreEqual(endDeflectionExpected, endDeflectionComputed));
        }
    }
}
