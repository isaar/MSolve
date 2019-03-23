using System;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Feti;
using ISAAR.MSolve.Solvers.Iterative;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: Ensure that the global dof ordering is the same as the one FETI uses
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Feti
{
    public class ExactPcpgResidualCalculator : IExactResidualCalculator
    {
        private readonly Func<Model_v2, ISolver_v2, IStaticProvider_v2> createProblemProvider;
        private readonly IDofOrderer originalDofOrderer;
        private readonly Model_v2 singleSubdomainModel;
        private readonly Subdomain_v2 singleSubdomain;
        private IMatrixView globalStiffness;
        private IVectorView globalForces;

        public ExactPcpgResidualCalculator(Model_v2 singleSubdomainModel, IDofOrderer originalDofOrderer, 
            Func<Model_v2, ISolver_v2, IStaticProvider_v2> createProblemProvider)
        {
            this.singleSubdomainModel = singleSubdomainModel;
            this.singleSubdomain = singleSubdomainModel.Subdomains.First();
            this.originalDofOrderer = originalDofOrderer;
            this.createProblemProvider = createProblemProvider;
        }

        public double CalculateExactResidualNorm(IVectorView globalDisplacements)
            => globalForces.Subtract(globalStiffness.Multiply(globalDisplacements)).Norm2();

        public void BuildLinearSystem()
        {
            // PcgSolver uses CSR matrices which are efficient for calculating f-K*u
            var pcgBuilder = new PcgAlgorithm.Builder();
            pcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(1); // No need to solve though.
            var solverBuilder = new PcgSolver.Builder();
            solverBuilder.DofOrderer = originalDofOrderer;
            solverBuilder.PcgAlgorithm = pcgBuilder.Build();
            PcgSolver solver = solverBuilder.BuildSolver(singleSubdomainModel);

            // Let MSolve follow the usual analysis routine, to create all necessary data structures. 
            IStaticProvider_v2 problemProvider = createProblemProvider(singleSubdomainModel, solver);
            var linearAnalyzer = new LinearAnalyzer_v2(singleSubdomainModel, solver, problemProvider);
            var staticAnalyzer = new StaticAnalyzer_v2(singleSubdomainModel, solver, problemProvider, linearAnalyzer);
            staticAnalyzer.Initialize();
            try
            {
                staticAnalyzer.Solve();
            }
            catch (IterativeSolverNotConvergedException)
            { }

            // Extract the global matrix and rhs
            ILinearSystem_v2 linearSystem = solver.LinearSystems.First().Value;
            globalStiffness = linearSystem.Matrix;
            globalForces = linearSystem.RhsVector;
        }
    }
}
