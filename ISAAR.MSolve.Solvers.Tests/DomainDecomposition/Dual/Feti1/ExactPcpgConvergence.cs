using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Projection;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;
using ISAAR.MSolve.Solvers.Iterative;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: Ensure that the global dof ordering is the same as the one FETI uses
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual
{
    public class ExactPcpgConvergence : IFetiPcgConvergence
    {
        private readonly Feti1Solver fetiSolver;
        private readonly IMatrixView globalStiffness;
        private readonly IVectorView globalForces;
        private readonly double globalForcesNorm;

        /// <summary>
        /// Could be null if another interface problem solver is used.
        /// </summary>
        private readonly Feti1ProjectedInterfaceProblemSolver interfaceProblemSolver; //TODO: avoid leaving this null.

        private readonly Feti1Projection projection; // only needed when separating the lagranges during PCG
        private readonly Vector lagrangesParticular; // only needed when separating the lagranges during PCG

        public ExactPcpgConvergence(IMatrixView globalStiffness, IVectorView globalForces, double globalForcesNorm, 
            Feti1Solver fetiSolver, Feti1ProjectedInterfaceProblemSolver interfaceProblemSolver,
            Feti1Projection projection, Vector lagrangesParticular)
        {
            this.globalStiffness = globalStiffness;
            this.globalForces = globalForces;
            this.globalForcesNorm = globalForcesNorm;
            this.fetiSolver = fetiSolver;
            this.interfaceProblemSolver = interfaceProblemSolver;
            this.projection = projection;
            this.lagrangesParticular = lagrangesParticular;
        }

        public double EstimateResidualNormRatio(PcgAlgorithmBase pcg)
        {
            var lagrangesBar = Vector.CreateZero(pcg.Solution.Length);
            lagrangesBar.CopyFrom(pcg.Solution);
            Vector lagranges = interfaceProblemSolver.CombineLagrangeMultipliers(lagrangesParticular, lagrangesBar, projection);
            return CalcExactResidualNorm(lagranges);
        }

        public double EstimateResidualNormRatio(IVectorView lagrangeMultipliers, IVectorView projectedPrecondResidual)
        {
            var lagrangesDense = Vector.CreateZero(lagrangeMultipliers.Length);
            lagrangesDense.CopyFrom(lagrangeMultipliers);
            return CalcExactResidualNorm(lagrangesDense);
        }

        public void Initialize(PcgAlgorithmBase pcg) { } // Do nothing

        //TODO: This is only used by the corresponding PCPG convergence criteria, which itself is for testing purposes only. 
        //      It should be moved to another class. If private members of this class are needed, then reflection should be used
        //      to access them.
        private double CalcExactResidualNorm(Vector lagranges)
        {
            Vector disconnectedDisplacements = fetiSolver.CalcDisconnectedDisplacements();
            Vector rbmCoeffs = fetiSolver.CalcRigidBodyModesCoefficients(disconnectedDisplacements, lagranges);
            var subdomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var idDisplacements in fetiSolver.CalcActualDisplacements(lagranges, rbmCoeffs))
            {
                subdomainDisplacements[idDisplacements.Key] = idDisplacements.Value;
            }
            Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements(subdomainDisplacements);
            return globalForces.Subtract(globalStiffness.Multiply(globalDisplacements)).Norm2() / globalForcesNorm;
        }

        public class Factory : IFetiPcgConvergenceFactory
        {
            private readonly Func<Model_v2, ISolver_v2, IStaticProvider_v2> createProblemProvider;
            private readonly double globalForcesNorm;
            private readonly IDofOrderer originalDofOrderer;
            private readonly Model_v2 singleSubdomainModel;
            private readonly Subdomain_v2 singleSubdomain;

            public Factory(Model_v2 singleSubdomainModel, IDofOrderer originalDofOrderer,
                Func<Model_v2, ISolver_v2, IStaticProvider_v2> createProblemProvider)
            {
                this.singleSubdomainModel = singleSubdomainModel;
                this.singleSubdomain = singleSubdomainModel.Subdomains.First();
                this.originalDofOrderer = originalDofOrderer;
                this.createProblemProvider = createProblemProvider;
            }

            public Feti1Solver FetiSolver { get; set; }

            /// <summary>
            /// Could be null if another interface problem solver is used.
            /// </summary>
            public Feti1ProjectedInterfaceProblemSolver InterfaceProblemSolver { get; set; }

            public IFetiPcgConvergence CreateConvergenceStrategy(double globalForcesNorm)
            {
                (IMatrixView globalStiffness, IVectorView globalForces) = BuildGlobalLinearSystem();

                if (FetiSolver == null) throw new InvalidOperationException(
                    "The FetiSolver property of this object must be set first");

                // Use reflection to access private members of FETI Solver
                var projection = (Feti1Projection)typeof(Feti1Solver).GetField("projection",
                    BindingFlags.NonPublic | BindingFlags.Instance).GetValue(FetiSolver);

                // λ0 = Q * G * inv(G^T * Q * G) * e
                Vector lagrangesParticular = projection.CalcParticularLagrangeMultipliers(FetiSolver.CalcRigidBodyModesWork());

                return new ExactPcpgConvergence(globalStiffness, globalForces, globalForcesNorm, 
                    FetiSolver, InterfaceProblemSolver, projection, lagrangesParticular);
            }

            //TODO: perhaps this should be done in Initialize(). Nope, Initialize is not called when using PCPG.
            private (IMatrixView globalStiffness, IVectorView globalForces) BuildGlobalLinearSystem()
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
                return (linearSystem.Matrix, linearSystem.RhsVector);
            }
        }
    }
}
