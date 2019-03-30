using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Preconditioning;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.Feti1
{
    /// <summary>
    /// Tests the correctness of the FETI-1 solver, with lumped preconditioner in a simple homogeneous problem.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SimplePlateTest
    {
        private const int singleSubdomainID = 0;

        [Fact]
        internal static void Run()
        {
            int numElementsX = 20, numElementsY = 10;
            double factorizationTolerance = 1E-4; // Defining the rigid body modes is very sensitive to this. TODO: The user shouldn't have to specify such a volatile parameter.
            double equalityTolerance = 1E-7;
            IVectorView expectedDisplacements = SolveModelWithoutSubdomains(numElementsX, numElementsY);
            (IVectorView computedDisplacements, DualSolverLogger logger) = 
                SolveModelWithSubdomains(numElementsX, numElementsY, factorizationTolerance);
            Debug.WriteLine($"Iterations: {logger.PcgIterations}");
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
            builder.NumTotalElementsX = numElementsX;
            builder.NumTotalElementsY = numElementsY;
            builder.YoungModulus = 2.1E7;
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerLeftCorner, DOFType.X, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerLeftCorner, DOFType.Y, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerRightCorner, DOFType.Y, 0.0);
            builder.PrescribeDistributedLoad(Uniform2DModelBuilder.BoundaryRegion.RightSide, DOFType.X, 100.0);

            return builder.BuildModel();
        }

        private static Model_v2 CreateSingleSubdomainModel(int numElementsX, int numElementsY)
        {
            // Replace the existing subdomains with a single one 
            Model_v2 model = CreateModel(numElementsX, numElementsY);
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain_v2(singleSubdomainID);
            model.SubdomainsDictionary.Add(singleSubdomainID, subdomain);
            foreach (Element_v2 element in model.Elements) subdomain.Elements.Add(element);
            return model;
        }

        private static (IVectorView U, DualSolverLogger logger) SolveModelWithSubdomains(int numElementsX, int numElementsY,
            double factorizationTolerance)
        {
            Model_v2 multiSubdomainModel = CreateModel(numElementsX, numElementsY);

            // Solver
            var solverBuilder = new Feti1Solver.Builder(factorizationTolerance);
            solverBuilder.PreconditionerFactory = new Feti1LumpedPreconditioner.Factory();
            solverBuilder.ProblemIsHomogeneous = true;
            Feti1Solver fetiSolver = solverBuilder.BuildSolver(multiSubdomainModel);

            // Linear static analysis
            var provider = new ProblemStructural_v2(multiSubdomainModel, fetiSolver);
            var childAnalyzer = new LinearAnalyzer_v2(multiSubdomainModel, fetiSolver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(multiSubdomainModel, fetiSolver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Gather the global displacements
            var sudomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var ls in fetiSolver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
            return (fetiSolver.GatherGlobalDisplacements(sudomainDisplacements), fetiSolver.Logger);
        }

        private static IVectorView SolveModelWithoutSubdomains(int numElementsX, int numElementsY)
        {
            Model_v2 model = CreateSingleSubdomainModel(numElementsX, numElementsY);

            // Solver
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);

            // Linear static analysis
            var provider = new ProblemStructural_v2(model, solver);
            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[singleSubdomainID].Solution;
        }
    }
}
