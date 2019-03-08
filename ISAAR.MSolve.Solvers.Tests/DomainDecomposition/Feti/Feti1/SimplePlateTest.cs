using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Feti;
using ISAAR.MSolve.Solvers.DomainDecomposition.Feti1;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Feti.Feti1
{
    /// <summary>
    /// Tests the correctness of the FETI-1 solver, with lumped preconditioner in a simple homogeneous problem.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SimplePlateTest
    {
        private const int subdomainIDsStart = 0;

        [Fact]
        internal static void Run()
        {
            int numElementsX = 20, numElementsY = 10;
            double factorizationTolerance = 1E-4; // Defining the rigid body modes is very sensitive to this. TODO: The user shouldn't have to specify such a volatile parameter.
            double equalityTolerance = 1E-7;
            var logger = new FetiLogger();
            IVectorView expectedDisplacements = SolveModelWithoutSubdomains(numElementsX, numElementsY);
            IVectorView computedDisplacements = SolveModelWithSubdomains(numElementsX, numElementsY, factorizationTolerance, 
                logger, true);
            Debug.WriteLine($"Iterations: {logger.PcpgIterations}");
            Assert.True(expectedDisplacements.Equals(computedDisplacements, equalityTolerance));
        }

        private static Model_v2 CreateModel(int numElementsX, int numElementsY)
        {
            // if numElementsX = numElementsY = 2:
            // 6 ----- 7 ----- 8  -->
            // |  (2)  |  (3)  |
            // |       |       |
            // 3 ----- 4 ----- 5  -->
            // |  (0)  |  (1)  |
            // |       |       |
            // 0 ----- 1 ----- 2  -->
            // Δ               Δ    
            // -               o

            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = 2.0;
            builder.DomainLengthY = 2.0;
            builder.NumSubdomainsX = 2;
            builder.NumSubdomainsY = 2;
            builder.NumElementsPerSubdomainX = numElementsX / builder.NumSubdomainsX;
            builder.NumElementsPerSubdomainY = numElementsY / builder.NumSubdomainsY;
            builder.YoungModulus = 2.1E7;
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerLeftCorner, DOFType.X, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerLeftCorner, DOFType.Y, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerRightCorner, DOFType.Y, 0.0);
            builder.PrescribeDistributedLoad(Uniform2DModelBuilder.BoundaryRegion.RightSide, DOFType.X, 100.0);

            return builder.CreateModel();
        }

        private static Model_v2 CreateSingleSubdomainModel(int numElementsX, int numElementsY)
        {
            // Replace the existing subdomains with a single one 
            Model_v2 model = CreateModel(numElementsX, numElementsY);
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain_v2(subdomainIDsStart);
            model.SubdomainsDictionary.Add(subdomainIDsStart, subdomain);
            foreach (Element_v2 element in model.Elements) subdomain.Elements.Add(element);
            return model;
        }

        private static IVectorView SolveModelWithSubdomains(int numElementsX, int numElementsY, double factorizationTolerance,
            FetiLogger logger, bool exactResidual)
        {
            // Solver
            Model_v2 multiSubdomainModel = CreateModel(numElementsX, numElementsY);
            var solverBuilder = new Feti1Solver.Builder(factorizationTolerance);
            solverBuilder.Logger = logger;

            // If PCPG needs to use the exact residual (only useful for testing) 
            if (exactResidual)
            {
                Model_v2 singleSubdomainModel = CreateSingleSubdomainModel(numElementsX, numElementsY);
                var exactResidualCalculator = new ExactPcpgResidualCalculator(singleSubdomainModel, 
                    solverBuilder.DofOrderer, (model, solver) => new ProblemStructural_v2(model, solver));
                exactResidualCalculator.BuildLinearSystem();
                solverBuilder.PcpgExactResidual = exactResidualCalculator;
            }

            // Structural problem provider
            Feti1Solver fetiSolver = solverBuilder.BuildSolver(multiSubdomainModel);
            var provider = new ProblemStructural_v2(multiSubdomainModel, fetiSolver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(multiSubdomainModel, fetiSolver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(multiSubdomainModel, fetiSolver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Gather the global displacements
            var sudomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var ls in fetiSolver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
            return fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);
        }

        private static IVectorView SolveModelWithoutSubdomains(int numElementsX, int numElementsY)
        {
            Model_v2 model = CreateSingleSubdomainModel(numElementsX, numElementsY);

            // Solver
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural_v2(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[subdomainIDsStart].Solution;
        }
    }
}
