using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.Feti1
{
    /// <summary>
    /// Tests from Papagiannakis bachelor thesis (NTUA 2011), p. 134 - 147
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class PapagiannakisFetiDPTests2D
    {
        public enum Precond { Dirichlet, DirichletDiagonal, Lumped }
        public enum Residual { Approximate, Exact }

        private const int singleSubdomainID = 0;
        private const int maxIterations = 1000;

        //TODO: Exact residual calculation is not implemented yet.
        [Theory]
        // Homogeneous problem
        [InlineData(1.0, Precond.Dirichlet, Residual.Exact, 11)]
        [InlineData(1.0, Precond.DirichletDiagonal, Residual.Exact, 14)]
        [InlineData(1.0, Precond.Lumped, Residual.Exact, 18)]

        // Stiffness ratio = 1E-2
        [InlineData(1E-2, Precond.Dirichlet, Residual.Exact, 12)]
        [InlineData(1E-2, Precond.DirichletDiagonal, Residual.Exact, 16)]
        [InlineData(1E-2, Precond.Lumped, Residual.Exact, 21)]


        // Stiffness ratio = 1E-3
        [InlineData(1E-3, Precond.Dirichlet, Residual.Exact, 13)]
        [InlineData(1E-3, Precond.DirichletDiagonal, Residual.Exact, 20)]
        [InlineData(1E-3, Precond.Lumped, Residual.Exact, 22)]

        // Stiffness ratio = 1E-4
        [InlineData(1E-4, Precond.Dirichlet, Residual.Exact, 14)]
        [InlineData(1E-4, Precond.DirichletDiagonal, Residual.Exact, 22)]
        [InlineData(1E-4, Precond.Lumped, Residual.Exact, 26)]

        // Stiffness ratio = 1E-5
        [InlineData(1E-5, Precond.Dirichlet, Residual.Exact, 14)]
        [InlineData(1E-5, Precond.DirichletDiagonal, Residual.Exact, 23)]
        [InlineData(1E-5, Precond.Lumped, Residual.Exact, 30)]

        // Stiffness ratio = 1E-6
        [InlineData(1E-6, Precond.Dirichlet, Residual.Exact, 15)]
        [InlineData(1E-6, Precond.DirichletDiagonal, Residual.Exact, 27)]
        [InlineData(1E-6, Precond.Lumped, Residual.Exact, 33)]
        public static void Run(double stiffnessRatio, Precond precond, Residual convergence, int iterExpected)
        {
            double pcpgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains(stiffnessRatio);
            (IVectorView ddDisplacements, DualSolverLogger logger, int numUniqueGlobalDofs, int numExtenedDomainDofs) =
                SolveModelWithSubdomains(stiffnessRatio, precond, convergence, pcpgConvergenceTol);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            Assert.Equal(882, numUniqueGlobalDofs);    // 882 includes constrained and free dofs
            Assert.Equal(1056, numExtenedDomainDofs); // 1056 includes constrained and free dofs
            Assert.Equal(190, logger.NumLagrangeMultipliers);

            // The error is provided in the reference solution the, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 5);

            // Allow a tolerance: It is ok if my solver is better or off by 1 iteration 
            Assert.InRange(logger.PcgIterations, 1, iterExpected + 1); // the upper bound is inclusive!
        }

        private static Model CreateModel(double stiffnessRatio)
        {
            // Subdomains:
            // /|
            // /||-------|-------|-------|-------|  
            // /||  (4)  |  (5)  |  (6)  |  (7)  |
            // /||   E1  |   E0  |   E0  |   E0  |
            // /||-------|-------|-------|-------|  
            // /||  (0)  |  (1)  |  (2)  |  (3)  |
            // /||   E1  |   E0  |   E0  |   E0  |
            // /||-------|-------|-------|-------|
            // /|

            double E0 = 2.1E7;
            double E1 = stiffnessRatio * E0;

            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = 3.0;
            builder.DomainLengthY = 1.5;
            builder.NumSubdomainsX = 4;
            builder.NumSubdomainsY = 2;
            builder.NumTotalElementsX = 20;
            builder.NumTotalElementsY = 20;
            builder.YoungModuliOfSubdomains = new double[,] { { E1, E0, E0, E0 }, { E1, E0, E0, E0 } };
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationY, 0.0);
            builder.DistributeLoadAtNodes(Uniform2DModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationY, 100.0);

            return builder.BuildModel();
        }

        private static Model CreateSingleSubdomainModel(double stiffnessRatio)
        {
            // Replace the existing subdomains with a single one 
            Model model = CreateModel(stiffnessRatio);
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain(singleSubdomainID);
            model.SubdomainsDictionary.Add(singleSubdomainID, subdomain);
            foreach (Element element in model.Elements) subdomain.Elements.Add(element);
            return model;
        }

        private static (IVectorView globalDisplacements, DualSolverLogger logger, int numUniqueGlobalDofs, int numExtenedDomainDofs)
            SolveModelWithSubdomains(double stiffnessRatio, Precond precond, Residual residualConvergence, 
                double pcgConvergenceTolerance)
        {
            // Model
            Model multiSubdomainModel = CreateModel(stiffnessRatio);

            // Corner nodes
            var cornerNodesOfEachSubdomain = new Dictionary<int, INode[]>();
            foreach (Subdomain subdomain in multiSubdomainModel.Subdomains)
            {
                subdomain.DefineNodesFromElements(); //TODO: This will also be called by the analyzer.
                int[] cornerNodeIDs = subdomain.GetCornerNodes();
                var cornerNodes = new List<INode>();
                foreach (int id in cornerNodeIDs)
                {
                    Node node = multiSubdomainModel.NodesDictionary[id];
                    if (node.Constraints.Count == 0) cornerNodes.Add(node);
                }
                cornerNodesOfEachSubdomain[subdomain.ID] = cornerNodes.ToArray();
            }

            // Solver
            var solverBuilder = new FetiDPSolver.Builder(cornerNodesOfEachSubdomain);

            // Homogeneous/heterogeneous problem
            solverBuilder.ProblemIsHomogeneous = stiffnessRatio == 1.0;

            // Preconditioner
            if (precond == Precond.Lumped) solverBuilder.PreconditionerFactory = new LumpedPreconditioner.Factory();
            else if (precond == Precond.DirichletDiagonal)
            {
                solverBuilder.PreconditionerFactory = new DiagonalDirichletPreconditioner.Factory();
            }
            else solverBuilder.PreconditionerFactory = new DirichletPreconditioner.Factory();

            //TODO: This needs to be implemented for FETI-DP
            //// PCG may need to use the exact residual for the comparison with the expected values
            //bool residualIsExact = residualConvergence == Residual.Exact;
            //ExactPcpgConvergence.Factory exactResidualConvergence = null;
            //if (residualIsExact)
            //{
            //    exactResidualConvergence = new ExactPcpgConvergence.Factory(
            //        CreateSingleSubdomainModel(stiffnessRatio), solverBuilder.DofOrderer,
            //            (model, solver) => new ProblemStructural(model, solver));
            //}
            //if (residualIsExact) exactResidualConvergence.InterfaceProblemSolver = interfaceProblemSolver;

            FetiDPSolver fetiSolver = solverBuilder.BuildSolver(multiSubdomainModel);
            //if (residualIsExact) exactResidualConvergence.FetiSolver = fetiSolver;

            // Structural problem provider
            var provider = new ProblemStructural(multiSubdomainModel, fetiSolver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(multiSubdomainModel, fetiSolver, provider);
            var parentAnalyzer = new StaticAnalyzer(multiSubdomainModel, fetiSolver, provider, childAnalyzer);

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

            return (globalDisplacements, fetiSolver.Logger, numUniqueGlobalDofs, numExtenedDomainDofs);
        }

        private static IVectorView SolveModelWithoutSubdomains(double stiffnessRatio)
        {
            Model model = CreateSingleSubdomainModel(stiffnessRatio);

            // Solver
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[singleSubdomainID].Solution;
        }
    }
}
