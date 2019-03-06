using System;
using System.Collections.Generic;
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
            IVectorView expectedDisplacements = SolveModelWithoutSubdomains(numElementsX, numElementsY);
            IVectorView computedDisplacements = SolveModelWithSubdomains(numElementsX, numElementsY, factorizationTolerance);
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

        private static IVectorView SolveModelWithSubdomains(int numElementsX, int numElementsY, double factorizationTolerance)
        {
            Model_v2 model = CreateModel(numElementsX, numElementsY);

            // Solver
            var solver = new Feti1Solver(model, factorizationTolerance);

            // Structural problem provider
            var provider = new ProblemStructural_v2(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.GatherGlobalDisplacements();
        }

        private static IVectorView SolveModelWithoutSubdomains(int numElementsX, int numElementsY)
        {
            Model_v2 model = CreateModel(numElementsX, numElementsY);

            // Replace the existing subdomains with a single one 
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain_v2(subdomainIDsStart);
            model.SubdomainsDictionary.Add(subdomainIDsStart, subdomain);
            foreach (Element_v2 element in model.Elements) subdomain.Elements.Add(element);

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
