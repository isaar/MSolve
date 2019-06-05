using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Rigid body modes do not have to be computed each time the stiffness matrix changes. 
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    public class FetiDPSolver : ISolver
    {
        internal const string name = "FETI-DP Solver"; // for error messages
        private readonly Dictionary<int, INode[]> cornerNodesOfSubdomains;
        private readonly ICrosspointStrategy crosspointStrategy = new FullyRedundantConstraints();
        private readonly IDofOrderer dofOrderer;
        private readonly IFetiDPInterfaceProblemSolver interfaceProblemSolver;
        private readonly Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers;
        private readonly Dictionary<int, IFetiSubdomainMatrixManager> matrixManagersGeneral; //TODO: redesign. They are the same as above, but Dictionary is not covariant
        private readonly Dictionary<int, ISingleSubdomainLinearSystem> linearSystems;
        private readonly IStructuralModel model;
        private readonly IFetiPreconditionerFactory preconditionerFactory;
        private readonly bool problemIsHomogeneous;

        //TODO: fix the mess of Dictionary<int, ISubdomain>, List<ISubdomain>, Dictionary<int, Subdomain>, List<Subdomain>
        //      The concrete are useful for the preprocessor mostly, while analyzers, solvers need the interface versions.
        //      Lists are better in analyzers and solvers. Dictionaries may be better in the preprocessors.
        private readonly Dictionary<int, ISubdomain> subdomains;

        private FetiDPDofSeparator dofSeparator;
        private bool factorizeInPlace = true;
        private FetiDPFlexibilityMatrix flexibility;
        private bool isStiffnessModified = true;
        private FetiDPLagrangeMultipliersEnumerator lagrangeEnumerator;
        private IFetiPreconditioner preconditioner;
        private IStiffnessDistribution stiffnessDistribution;
        private FetiDPSubdomainGlobalMapping subdomainGlobalMapping;

        private FetiDPSolver(IStructuralModel model, Dictionary<int, INode[]> cornerNodesOfSubdomains,
            IFetiDPSubdomainMatrixManagerFactory matrixManagerFactory, IDofOrderer dofOrderer, 
            IFetiPreconditionerFactory preconditionerFactory, IFetiDPInterfaceProblemSolver interfaceProblemSolver, 
            bool problemIsHomogeneous)
        {
            // Model
            if (model.Subdomains.Count == 1) throw new InvalidSolverException(
                $"{name} cannot be used if there is only 1 subdomain");
            this.model = model;
            this.cornerNodesOfSubdomains = cornerNodesOfSubdomains;

            // Subdomains
            subdomains = new Dictionary<int, ISubdomain>();
            foreach (ISubdomain subdomain in model.Subdomains) subdomains[subdomain.ID] = subdomain;

            // Matrix managers and linear systems
            matrixManagers = new Dictionary<int, IFetiDPSubdomainMatrixManager>();
            matrixManagersGeneral = new Dictionary<int, IFetiSubdomainMatrixManager>();
            this.linearSystems = new Dictionary<int, ISingleSubdomainLinearSystem>();
            var externalLinearSystems = new Dictionary<int, ILinearSystem>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                var matrixManager = matrixManagerFactory.CreateMatricesManager(subdomain);
                matrixManagers[id] = matrixManager;
                matrixManagersGeneral[id] = matrixManager;
                this.linearSystems[id] = matrixManager.LinearSystem;
                externalLinearSystems[id] = matrixManager.LinearSystem;
            }
            LinearSystems = externalLinearSystems;

            this.dofOrderer = dofOrderer;
            this.preconditionerFactory = preconditionerFactory;

            // Interface problem
            this.interfaceProblemSolver = interfaceProblemSolver;

            // Homogeneous/heterogeneous problems
            this.problemIsHomogeneous = problemIsHomogeneous;
        }

        public IReadOnlyDictionary<int, ILinearSystem> LinearSystems { get; }
        public DualSolverLogger Logger { get; } = new DualSolverLogger();

        public Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider elementMatrixProvider)
        {
            var matrices = new Dictionary<int, IMatrix>();
            var matricesReadonly = new Dictionary<int, IMatrixView>();
            foreach (ISubdomain subdomain in model.Subdomains) //TODO: this must be done in parallel
            {
                IMatrix stiffness = matrixManagers[subdomain.ID].BuildGlobalMatrix(subdomain.FreeDofOrdering,
                    subdomain.Elements, elementMatrixProvider);
                matricesReadonly[subdomain.ID] = stiffness;
                matrices[subdomain.ID] = stiffness;
            }

            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            DetermineStiffnessDistribution(matricesReadonly);

            return matrices;
        }

        public Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr)> BuildGlobalSubmatrices(IElementMatrixProvider elementMatrixProvider)
        {
            var matrices = new Dictionary<int, (IMatrix Aff, IMatrixView Afc, IMatrixView Acf, IMatrixView Acc)>();
            var matricesReadonly = new Dictionary<int, IMatrixView>();
            foreach (ISubdomain subdomain in model.Subdomains) //TODO: this must be done in parallel
            {
                if (subdomain.ConstrainedDofOrdering == null)
                {
                    throw new InvalidOperationException("In order to build the matrices corresponding to constrained dofs of,"
                        + $" subdomain {subdomain.ID}, they must have been ordered first.");
                }
                (IMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) =
                    matrixManagers[subdomain.ID].BuildGlobalSubmatrices(subdomain.FreeDofOrdering,
                    subdomain.ConstrainedDofOrdering, subdomain.Elements, elementMatrixProvider);
                matrices[subdomain.ID] = (Kff, Kfc, Kcf, Kcc);
                matricesReadonly[subdomain.ID] = Kff;
            }

            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            DetermineStiffnessDistribution(matricesReadonly);

            return matrices;
        }

        //TODO: this and the fields should be handled by a class that handles dof mappings.
        public Dictionary<int, SparseVector> DistributeNodalLoads(Table<INode, IDofType, double> globalNodalLoads)
            => subdomainGlobalMapping.DistributeNodalLoads(subdomains, globalNodalLoads);

        //TODO: this and the fields should be handled by a class that handles dof mappings.
        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainDisplacements)
            => subdomainGlobalMapping.GatherGlobalDisplacements(subdomainDisplacements);

        public void HandleMatrixWillBeSet()
        {
            isStiffnessModified = true;
            foreach (IFetiDPSubdomainMatrixManager matrixManager in matrixManagers.Values) matrixManager.Clear();
            flexibility = null;
            preconditioner = null;
            interfaceProblemSolver.ClearCoarseProblemMatrix();

            //stiffnessDistribution = null; //WARNING: do not dispose of this. It is updated when BuildGlobalMatrix() is called.
        }

        public void Initialize()
        { }

        public Dictionary<int, Matrix> InverseSystemMatrixTimesOtherMatrix(Dictionary<int, IMatrixView> otherMatrix)
        {
            throw new NotImplementedException();
        }

        public void OrderDofs(bool alsoOrderConstrainedDofs)
        {
            // Order dofs
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderFreeDofs(model);
            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                matrixManagers[subdomain.ID].HandleDofOrderingWillBeModified();
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
                if (alsoOrderConstrainedDofs) subdomain.ConstrainedDofOrdering = dofOrderer.OrderConstrainedDofs(subdomain);

                // The next must done by the analyzer, so that subdomain.Forces is retained when doing back to back analyses.
                //subdomain.Forces = linearSystem.CreateZeroVector();
            }

            // Define boundary / internal dofs
            dofSeparator = new FetiDPDofSeparator();
            dofSeparator.SeparateDofs(model, cornerNodesOfSubdomains);
            dofSeparator.DefineCornerMappingMatrices(model, cornerNodesOfSubdomains);

            // Define lagrange multipliers and boolean matrices
            this.lagrangeEnumerator = new FetiDPLagrangeMultipliersEnumerator(crosspointStrategy, dofSeparator);
            if (problemIsHomogeneous) lagrangeEnumerator.DefineBooleanMatrices(model); // optimization in this case
            else lagrangeEnumerator.DefineLagrangesAndBooleanMatrices(model);

            // Log dof statistics
            if (Logger != null)
            {
                Logger.NumUniqueGlobalFreeDofs = model.GlobalDofOrdering.NumGlobalFreeDofs;
                Logger.NumExpandedDomainFreeDofs = 0;
                foreach (var subdomain in model.Subdomains)
                {
                    Logger.NumExpandedDomainFreeDofs += subdomain.FreeDofOrdering.NumFreeDofs;
                }
                Logger.NumLagrangeMultipliers = lagrangeEnumerator.NumLagrangeMultipliers;
                Logger.NumCornerDofs = dofSeparator.NumGlobalCornerDofs;
            }
        }

        public void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        public void Solve()
        {
            foreach (var linearSystem in linearSystems.Values)
            {
                if (linearSystem.SolutionConcrete == null) linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
            }

            // Separate the force vector
            var fr = new Dictionary<int, Vector>();
            var fbc = new Dictionary<int, Vector>();
            foreach (int s in subdomains.Keys)
            {
                int[] remainderDofs = dofSeparator.RemainderDofIndices[s];
                int[] cornerDofs = dofSeparator.CornerDofIndices[s];
                Vector f = linearSystems[s].RhsConcrete;
                fr[s] = f.GetSubvector(remainderDofs);
                fbc[s] = f.GetSubvector(cornerDofs);
            }

            if (isStiffnessModified)
            {
                // Separate the stiffness matrix
                foreach (int s in subdomains.Keys)
                {
                    IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];
                    int[] remainderDofs = dofSeparator.RemainderDofIndices[s];
                    int[] cornerDofs = dofSeparator.CornerDofIndices[s];
                    matrices.ExtractKrr(remainderDofs);
                    matrices.ExtractKcrKrc(cornerDofs, remainderDofs);
                    matrices.ExtractKcc(cornerDofs);
                }

                // Calculate the preconditioner before factorizing each subdomain's Kff 
                preconditioner = preconditionerFactory.CreatePreconditioner(stiffnessDistribution, dofSeparator,
                    lagrangeEnumerator, matrixManagersGeneral);

                // Factorize each subdomain's Krr
                foreach (int s in subdomains.Keys) matrixManagers[s].InvertKrr(true);

                // Define FETI-DP flexibility matrices
                flexibility = new FetiDPFlexibilityMatrix(dofSeparator, lagrangeEnumerator, matrixManagers);

                // Static condensation of remainder dofs (Schur complement).
                interfaceProblemSolver.CreateCoarseProblemMatrix(dofSeparator, matrixManagers);

                isStiffnessModified = false;
            }

            // Static condensation for the force vectors
            Vector globalFcStar = interfaceProblemSolver.CreateCoarseProblemRhs(dofSeparator, matrixManagers, fr, fbc);

            // Calculate the rhs vectors of the interface system
            Vector dr = CalcDisconnectedDisplacements(fr);
            double globalForcesNorm = CalcGlobalForcesNorm();

            // Solve the interface problem
            (Vector lagranges, Vector uc) = interfaceProblemSolver.SolveInterfaceProblem(flexibility, preconditioner, 
                globalFcStar, dr, globalForcesNorm, Logger);

            // Calculate the displacements of each subdomain
            Dictionary<int, Vector> actualDisplacements = CalcActualDisplacements(lagranges, uc,  fr);
            foreach (var idSystem in linearSystems) idSystem.Value.SolutionConcrete = actualDisplacements[idSystem.Key];
        }

        /// <summary>
        /// Does not mutate this object.
        /// </summary>
        internal Dictionary<int, Vector> CalcActualDisplacements(Vector lagranges, Vector cornerDisplacements, 
            Dictionary<int, Vector> fr)
        {
            var freeDisplacements = new Dictionary<int, Vector>();
            foreach (int s in subdomains.Keys)
            {
                IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];

                // ur[s] = inv(Krr[s]) * (fr[s] - Br[s]^T * lagranges - Krc[s] * Lc[s] * uc)
                Vector BrLambda = lagrangeEnumerator.BooleanMatrices[s].Multiply(lagranges, true);
                Vector KrcLcUc = dofSeparator.CornerBooleanMatrices[s].Multiply(cornerDisplacements);
                KrcLcUc = matrices.MultiplyKrcTimes(KrcLcUc);
                Vector temp = fr[s].Copy();
                temp.SubtractIntoThis(BrLambda);
                temp.SubtractIntoThis(KrcLcUc);
                Vector ur = matrices.MultiplyInverseKrrTimes(temp);

                // uf[s] = union(ur[s], ubc[s])
                // Remainder dofs
                var uf = Vector.CreateZero(subdomains[s].FreeDofOrdering.NumFreeDofs);
                int[] remainderDofs = dofSeparator.RemainderDofIndices[s];
                uf.CopyNonContiguouslyFrom(remainderDofs, ur);

                // Corner dofs: ubc[s] = Bc[s] * uc
                Vector ubc = dofSeparator.CornerBooleanMatrices[s].Multiply(cornerDisplacements);
                int[] cornerDofs = dofSeparator.CornerDofIndices[s];
                uf.CopyNonContiguouslyFrom(cornerDofs, ubc);

                freeDisplacements[s] = uf;
            }
            return freeDisplacements;
        }

        /// <summary>
        /// d = sum(Bs * generalInverse(Ks) * fs), where fs are the nodal forces applied to the dofs of subdomain s.
        /// Does not mutate this object.
        /// </summary>
        internal Vector CalcDisconnectedDisplacements(Dictionary<int, Vector> fr)
        {
            // dr = sum_over_s( Br[s] * inv(Krr[s]) * fr[s])
            var dr = Vector.CreateZero(lagrangeEnumerator.NumLagrangeMultipliers);
            foreach (int s in linearSystems.Keys)
            {
                SignedBooleanMatrixColMajor Br = lagrangeEnumerator.BooleanMatrices[s];
                Vector temp = matrixManagers[s].MultiplyInverseKrrTimes(fr[s]);
                temp = Br.Multiply(temp);
                dr.AddIntoThis(temp);
            }
            return dr;
        }

        private void BuildPreconditioner(Dictionary<int, Matrix> matricesKrr)
        {
            // Create the preconditioner. 
            //TODO: this should be done simultaneously with the factorizations to avoid duplicate factorizations.
            var stiffnessMatrices = new Dictionary<int, IMatrixView>();
            foreach (var idKrr in matricesKrr) stiffnessMatrices.Add(idKrr.Key, idKrr.Value);
            preconditioner = preconditionerFactory.CreatePreconditioner(stiffnessDistribution, dofSeparator,
                lagrangeEnumerator, null /*stiffnessMatrices*/);
        }

        /// <summary>
        /// Calculate the norm of the forces vector |f| = |K*u|. It is needed to check the convergence of PCG/PCPG.
        /// </summary>
        private double CalcGlobalForcesNorm()
        {
            //TODO: It would be better to do that using the global vector to avoid the homogeneous/heterogeneous averaging
            //      That would require the analyzer to build the global vector too though. Caution: we cannot take just 
            //      the nodal loads from the model, since the effect of other loads is only taken into account int 
            //      linearSystem.Rhs. Even if we could easily create the global forces vector, it might be wrong since 
            //      the analyzer may modify some of these loads, depending on time step, loading step, etc.
            var subdomainForces = new Dictionary<int, IVectorView>();
            foreach (var linearSystem in linearSystems.Values)
            {
                subdomainForces[linearSystem.Subdomain.ID] = linearSystem.RhsConcrete;
            }
            return subdomainGlobalMapping.CalculateGlobalForcesNorm(subdomainForces);
        }

        private void DetermineStiffnessDistribution(Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            if (problemIsHomogeneous)
            {
                stiffnessDistribution = new HomogeneousStiffnessDistribution(model, dofSeparator);
            }
            else
            {
                Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses = 
                    BoundaryDofLumpedStiffness.ExtractBoundaryDofLumpedStiffnesses(
                        dofSeparator.GlobalBoundaryDofs, stiffnessMatrices);
                stiffnessDistribution = new HeterogeneousStiffnessDistribution(model, dofSeparator, boundaryDofStiffnesses);
            }
            subdomainGlobalMapping = new FetiDPSubdomainGlobalMapping(model, dofSeparator, stiffnessDistribution);
        }

        public class Builder
        {
            private Dictionary<int, INode[]> cornerNodesOfEachSubdomain; //TODO: These should probably be a HashSet instead of array.
            private readonly IFetiDPSubdomainMatrixManagerFactory matrixManagerFactory;

            public Builder(Dictionary<int, INode[]> cornerNodesOfEachSubdomain, 
                IFetiDPSubdomainMatrixManagerFactory matrixManagerFactory)
            {
                this.cornerNodesOfEachSubdomain = cornerNodesOfEachSubdomain;
                this.matrixManagerFactory = matrixManagerFactory;
            }

            //TODO: We need to specify the ordering for remainder and possibly internal dofs, while IDofOrderer only works for free dofs.
            public IDofOrderer DofOrderer { get; set; } =
                new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public IFetiDPInterfaceProblemSolver InterfaceProblemSolver { get; set; } =
                new FetiDPInterfaceProblemSolver.Builder().Build();

            public IFetiPreconditionerFactory PreconditionerFactory { get; set; } = new LumpedPreconditioner.Factory();
            public bool ProblemIsHomogeneous { get; set; } = true;

            public FetiDPSolver BuildSolver(IStructuralModel model)
                => new FetiDPSolver(model, cornerNodesOfEachSubdomain, matrixManagerFactory, DofOrderer, PreconditionerFactory, 
                     InterfaceProblemSolver, ProblemIsHomogeneous);
        }
    }
}
