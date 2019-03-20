using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition;
using ISAAR.MSolve.Solvers.DomainDecomposition.Feti;
using ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Feti1;
using Xunit;

//TODO: the preconditioner could be provided as a parameter in [Theory] methods. At least the common code should be extracted.
//TODO: there is also data about the case without preconditioner and the case where Q=I, even in heterogeneous problems.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Feti.Feti1
{
    /// <summary>
    /// Tests from Papagiannakis bachelor thesis (NTUA 2011), p. 97 - 100
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class HomogeneousCantileverTests
    {
        private const int singleSubdomainID = 0;

        [Fact]
        public static void TestDiagonalDirichletPreconditioner()
        {
            double factorizationTol = 1E-3, pcpgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains();
            (IVectorView ddDisplacements, FetiLogger logger, int numUniqueGlobalDofs, int numExtenedDomainDofs) =
                SolveModelWithSubdomains(new Feti1DiagonalDirichletPreconditioner.Factory(), factorizationTol, 
                pcpgConvergenceTol, true);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            Assert.Equal(882, numUniqueGlobalDofs);    // 882 includes constrained and free dofs
            Assert.Equal(1056, numExtenedDomainDofs); // 1056 includes constrained and free dofs
            Assert.Equal(190, logger.NumLagrangeMultipliers);
            Assert.Equal(11, logger.PcpgIterations);

            // In the reference solution the error is 1.23E-9, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 8);
        }

        [Fact]
        public static void TestDirichletPreconditioner()
        {
            double factorizationTol = 1E-3, pcpgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains();
            (IVectorView ddDisplacements, FetiLogger logger, int numUniqueGlobalDofs, int numExtenedDomainDofs) =
                SolveModelWithSubdomains(new Feti1DirichletPreconditioner.Factory(), factorizationTol, pcpgConvergenceTol, true);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            Assert.Equal(882, numUniqueGlobalDofs);    // 882 includes constrained and free dofs
            Assert.Equal(1056, numExtenedDomainDofs); // 1056 includes constrained and free dofs
            Assert.Equal(190, logger.NumLagrangeMultipliers);
            Assert.Equal(10, logger.PcpgIterations);

            // In the reference solution the error is 2.59E-10, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 8);
        }

        [Fact]
        public static void TestLumpedPreconditioner()
        {
            double factorizationTol = 1E-3, pcpgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains();
            (IVectorView ddDisplacements, FetiLogger logger, int numUniqueGlobalDofs, int numExtenedDomainDofs) = 
                SolveModelWithSubdomains(new Feti1LumpedPreconditioner.Factory(), factorizationTol, pcpgConvergenceTol, true);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            Assert.Equal(882, numUniqueGlobalDofs);    // 882 includes constrained and free dofs
            Assert.Equal(1056, numExtenedDomainDofs); // 1056 includes constrained and free dofs
            Assert.Equal(190, logger.NumLagrangeMultipliers);

            // The reference solution needed 13 iterations but that depends on how ||f-A*0|| is calculated
            Assert.Equal(12, logger.PcpgIterations); 

            // In the reference solution the error is 5.5E-10, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 7);
        }

        private static Model_v2 CreateModel()
        {
            // Subdomains:
            // /|
            // /||-------|-------|-------|-------|  
            // /||  (4)  |  (5)  |  (6)  |  (7)  |
            // /||       |       |       |       |
            // /||-------|-------|-------|-------|  
            // /||  (0)  |  (1)  |  (2)  |  (3)  |
            // /||       |       |       |       |
            // /||-------|-------|-------|-------|
            // /|

            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = 3.0;
            builder.DomainLengthY = 1.5;
            builder.NumSubdomainsX = 4;
            builder.NumSubdomainsY = 2;
            builder.NumTotalElementsX = 20;
            builder.NumTotalElementsY = 20;
            builder.YoungModulus = 2.1E7;
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, DOFType.X, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, DOFType.Y, 0.0);
            builder.PrescribeDistributedLoad(Uniform2DModelBuilder.BoundaryRegion.RightSide, DOFType.Y, 100.0);

            return builder.BuildModel();
        }

        private static Model_v2 CreateSingleSubdomainModel()
        {
            // Replace the existing subdomains with a single one 
            Model_v2 model = CreateModel();
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain_v2(singleSubdomainID);
            model.SubdomainsDictionary.Add(singleSubdomainID, subdomain);
            foreach (Element_v2 element in model.Elements) subdomain.Elements.Add(element);
            return model;
        }

        private static (IVectorView globalDisplacements, FetiLogger logger, int numUniqueGlobalDofs, int numExtenedDomainDofs) 
            SolveModelWithSubdomains(IFetiPreconditionerFactory preconditioning, double factorizationTolerance, 
                double pcpgConvergenceTolerance, bool exactResidual)
        {
            // Model
            Model_v2 multiSubdomainModel = CreateModel();

            // Solver
            var solverBuilder = new Feti1Solver.Builder(factorizationTolerance);
            solverBuilder.IsProblemHomogeneous = true;
            solverBuilder.PcpgConvergenceTolerance = pcpgConvergenceTolerance;
            solverBuilder.PreconditionerFactory = preconditioning;
            solverBuilder.Logger = new FetiLogger();

            // PCPG needs to use the exact residual for the comparison with the expected values
            if (exactResidual)
            {
                var exactResidualCalculator = new ExactPcpgResidualCalculator(CreateSingleSubdomainModel(),
                   solverBuilder.DofOrderer, (model, solver) => new ProblemStructural_v2(model, solver));
                exactResidualCalculator.BuildLinearSystem();
                solverBuilder.PcpgExactResidual = exactResidualCalculator;
            }
            Feti1Solver fetiSolver = solverBuilder.BuildSolver(multiSubdomainModel);

            // Structural problem provider
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
            Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);

            // Other stats
            int numUniqueGlobalDofs = multiSubdomainModel.Nodes.Count * 2;
            int numExtenedDomainDofs = 0;
            foreach (var subdomain in multiSubdomainModel.Subdomains) numExtenedDomainDofs += subdomain.Nodes.Count * 2;

            return (globalDisplacements, solverBuilder.Logger, numUniqueGlobalDofs, numExtenedDomainDofs);
        }

        private static IVectorView SolveModelWithoutSubdomains()
        {
            Model_v2 model = CreateSingleSubdomainModel();

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

            return solver.LinearSystems[singleSubdomainID].Solution;
        }
    }
}
