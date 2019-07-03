using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP
{
    /// <summary>
    /// Tests from Papagiannakis bachelor thesis (NTUA 2011), p. 134 - 147
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class PapagiannakisFetiDPTests2D
    {
        public enum Precond { Dirichlet, DirichletDiagonal, Lumped }
        public enum Residual { Approximate, Exact }

        private const double domainLengthX = 3.0, domainLengthY = 1.5;
        private const int singleSubdomainID = 0;
        private const int maxIterations = 1000;

        //TODO: Exact residual calculation is not implemented yet. Therefore the expected iterations cannot be tested 
        //      accurately using the approximate residual yet.
        [Theory]
        // Homogeneous problem
        [InlineData(1.0, Precond.Dirichlet, Residual.Approximate, 11)]
        [InlineData(1.0, Precond.DirichletDiagonal, Residual.Approximate, 14)]
        [InlineData(1.0, Precond.Lumped, Residual.Approximate, 18)]

        // Stiffness ratio = 1E-2
        [InlineData(1E-2, Precond.Dirichlet, Residual.Approximate, 12)]
        [InlineData(1E-2, Precond.DirichletDiagonal, Residual.Approximate, 16)]
        [InlineData(1E-2, Precond.Lumped, Residual.Approximate, 21)]

        // Stiffness ratio = 1E-3
        [InlineData(1E-3, Precond.Dirichlet, Residual.Approximate, 13)]
        [InlineData(1E-3, Precond.DirichletDiagonal, Residual.Approximate, 20)]
        [InlineData(1E-3, Precond.Lumped, Residual.Approximate, 22)]

        // Stiffness ratio = 1E-4
        [InlineData(1E-4, Precond.Dirichlet, Residual.Approximate, 14)]
        [InlineData(1E-4, Precond.DirichletDiagonal, Residual.Approximate, 22)]
        [InlineData(1E-4, Precond.Lumped, Residual.Approximate, 26)]

        // Stiffness ratio = 1E-5
        [InlineData(1E-5, Precond.Dirichlet, Residual.Approximate, 14)]
        [InlineData(1E-5, Precond.DirichletDiagonal, Residual.Approximate, 23)]
        [InlineData(1E-5, Precond.Lumped, Residual.Approximate, 30)]

        // Stiffness ratio = 1E-6
        [InlineData(1E-6, Precond.Dirichlet, Residual.Approximate, 15)]
        [InlineData(1E-6, Precond.DirichletDiagonal, Residual.Approximate, 27)]
        [InlineData(1E-6, Precond.Lumped, Residual.Approximate, 33)]
        public static void Run(double stiffnessRatio, Precond precond, Residual convergence, int iterExpected)
        {
            double pcgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains(stiffnessRatio);
            (IVectorView ddDisplacements, SolverLogger logger) =
                SolveModelWithSubdomains(stiffnessRatio, precond, convergence, pcgConvergenceTol);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            int analysisStep = 0;
            Assert.Equal(140, logger.GetNumDofs(analysisStep, "Lagrange multipliers"));
            Assert.Equal(20, logger.GetNumDofs(analysisStep, "Corner dofs"));

            // The error is provided in the reference solution the, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 6);

            // Allow some tolerance for the iterations:
            int maxIterationsForApproximateResidual = (int)Math.Ceiling(1.0 * iterExpected);
            int pcgIterations = logger.GetNumIterationsOfIterativeAlgorithm(analysisStep);
            Assert.InRange(pcgIterations, 1, maxIterationsForApproximateResidual); // the upper bound is inclusive!
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
            builder.DomainLengthX = domainLengthX;
            builder.DomainLengthY = domainLengthY;
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

        private static (IVectorView globalDisplacements, SolverLogger logger) SolveModelWithSubdomains(double stiffnessRatio,
            Precond precond, Residual residualConvergence, double pcgConvergenceTolerance)
        {
            // Model
            Model multiSubdomainModel = CreateModel(stiffnessRatio);

            // Corner nodes
            double meshTol = 1E-6;
            var cornerNodesOfEachSubdomain = new Dictionary<int, HashSet<INode>>();
            foreach (Subdomain subdomain in multiSubdomainModel.Subdomains)
            {
                subdomain.DefineNodesFromElements(); //TODO: This will also be called by the analyzer.
                INode[] corners = CornerNodeUtilities.FindCornersOfRectangle2D(subdomain);
                var cornerNodes = new HashSet<INode>();
                foreach (INode node in corners)
                {
                    if (node.Constraints.Count > 0) continue;
                    if ((Math.Abs(node.X - domainLengthX) <= meshTol) && (Math.Abs(node.Y) <= meshTol)) continue;
                    if ((Math.Abs(node.X - domainLengthX) <= meshTol) && (Math.Abs(node.Y - domainLengthY) <= meshTol)) continue;
                    cornerNodes.Add(node);
                }
                cornerNodesOfEachSubdomain[subdomain.ID] = cornerNodes;
            }

            // Solver
            var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory(new OrderingAmdSuiteSparse());
            //var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
            //var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();
            var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodesOfEachSubdomain);
            var solverBuilder = new FetiDPSolver.Builder(cornerNodeSelection, fetiMatrices);
            solverBuilder.ProblemIsHomogeneous = stiffnessRatio == 1.0;
            //solverBuilder.ProblemIsHomogeneous = false;

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

            // Specify solver for the interface problem
            var interfaceSolverBuilder = new FetiDPInterfaceProblemSolver.Builder();
            interfaceSolverBuilder.PcgConvergenceStrategyFactory = new ApproximateResidualConvergence.Factory();
            interfaceSolverBuilder.PcgConvergenceTolerance = pcgConvergenceTolerance;
            solverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();

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

            return (globalDisplacements, fetiSolver.Logger);
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
