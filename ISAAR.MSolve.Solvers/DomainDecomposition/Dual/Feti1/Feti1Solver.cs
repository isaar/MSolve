using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Projection;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Rigid body modes do not have to be computed each time the stiffness matrix changes. 
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1
{
    public class Feti1Solver : ISolver
    {
        internal const string name = "FETI Level-1 Solver"; // for error messages
        private readonly ICrosspointStrategy crosspointStrategy = new FullyRedundantConstraints();
        private readonly IDofOrderer dofOrderer;
        private readonly Dictionary<int, double> factorPivotTolerances;
        private readonly IFeti1InterfaceProblemSolver interfaceProblemSolver;
        private readonly Dictionary<int, ISingleSubdomainLinearSystem> linearSystems;
        private readonly Dictionary<int, IFeti1SubdomainMatrixManager> matrixManagers;
        private readonly Dictionary<int, IFetiSubdomainMatrixManager> matrixManagersGeneral; //TODO: redesign. They are the same as above, but Dictionary is not covariant
        private readonly IStructuralModel model;
        //private readonly PdeOrder pde; // Instead the user explicitly sets Q.
        private readonly IFetiPreconditionerFactory preconditionerFactory;
        private readonly bool problemIsHomogeneous;
        private readonly bool projectionMatrixQIsIdentity;

        //TODO: fix the mess of Dictionary<int, ISubdomain>, List<ISubdomain>, Dictionary<int, Subdomain>, List<Subdomain>
        //      The concrete are useful for the preprocessor mostly, while analyzers, solvers need the interface versions.
        //      Lists are better in analyzers and solvers. Dictionaries may be better in the preprocessors.
        private readonly Dictionary<int, ISubdomain> subdomains;

        private Feti1DofSeparator dofSeparator;
        private bool factorizeInPlace = true;
        private Feti1FlexibilityMatrix flexibility;
        private bool isStiffnessModified = true;
        private Feti1LagrangeMultipliersEnumerator lagrangeEnumerator;
        private IFetiPreconditioner preconditioner;
        private Feti1Projection projection;
        private IStiffnessDistribution stiffnessDistribution;
        private Feti1SubdomainGlobalMapping subdomainGlobalMapping;

        private Feti1Solver(IStructuralModel model, IFeti1SubdomainMatrixManagerFactory matrixManagerFactory, 
            IDofOrderer dofOrderer, Dictionary<int, double> factorPivotTolerances, 
            IFetiPreconditionerFactory preconditionerFactory, IFeti1InterfaceProblemSolver interfaceProblemSolver, 
            bool problemIsHomogeneous, bool projectionMatrixQIsIdentity)
        {
            // Model
            if (model.Subdomains.Count == 1) throw new InvalidSolverException(
                $"{name} cannot be used if there is only 1 subdomain");
            this.model = model;

            // Subdomains
            subdomains = new Dictionary<int, ISubdomain>();
            foreach (ISubdomain subdomain in model.Subdomains) subdomains[subdomain.ID] = subdomain;

            // Matrix managers and linear systems
            matrixManagers = new Dictionary<int, IFeti1SubdomainMatrixManager>();
            matrixManagersGeneral = new Dictionary<int, IFetiSubdomainMatrixManager>();
            this.linearSystems = new Dictionary<int, ISingleSubdomainLinearSystem>();
            var externalLinearSystems = new Dictionary<int, ILinearSystem>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int s = subdomain.ID;
                var matrixManager = matrixManagerFactory.CreateMatricesManager(subdomain);
                matrixManagers[s] = matrixManager;
                matrixManagersGeneral[s] = matrixManager;
                this.linearSystems[s] = matrixManager.LinearSystem;
                externalLinearSystems[s] = matrixManager.LinearSystem;
            }
            LinearSystems = externalLinearSystems;

            this.dofOrderer = dofOrderer;
            this.factorPivotTolerances = factorPivotTolerances;
            this.preconditionerFactory = preconditionerFactory;

            // Interface problem
            this.interfaceProblemSolver = interfaceProblemSolver;

            // Homogeneous/heterogeneous problems
            this.problemIsHomogeneous = problemIsHomogeneous;
            this.projectionMatrixQIsIdentity = projectionMatrixQIsIdentity;
        }

        public IReadOnlyDictionary<int, ILinearSystem> LinearSystems { get; }
        public SolverLogger Logger { get; } = new SolverLogger(name);
        public string Name => name;

        public Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider elementMatrixProvider)
        {
            var watch = new Stopwatch();
            watch.Start();
            var matrices = new Dictionary<int, IMatrix>();
            var matricesReadonly = new Dictionary<int, IMatrixView>();
            foreach (ISubdomain subdomain in model.Subdomains) //TODO: this must be done in parallel
            {
                IMatrix stiffness = matrixManagers[subdomain.ID].BuildGlobalMatrix(subdomain.FreeDofOrdering,
                    subdomain.Elements, elementMatrixProvider);
                matricesReadonly[subdomain.ID] = stiffness;
                matrices[subdomain.ID] = stiffness;
            }
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);

            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            DetermineStiffnessDistribution(matricesReadonly);

            return matrices;
        }

        public Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr)> BuildGlobalSubmatrices(IElementMatrixProvider elementMatrixProvider)
        {
            var watch = new Stopwatch();
            watch.Start();
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
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);

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
            foreach (IFeti1SubdomainMatrixManager matrixManager in matrixManagers.Values) matrixManager.Clear();
            flexibility = null;
            preconditioner = null;
            projection = null;
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
            var watch = new Stopwatch();
            watch.Start();

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
            dofSeparator = new Feti1DofSeparator();
            dofSeparator.SeparateDofs(model);

            // Define lagrange multipliers and boolean matrices
            this.lagrangeEnumerator = new Feti1LagrangeMultipliersEnumerator(crosspointStrategy, dofSeparator);
            if (problemIsHomogeneous) lagrangeEnumerator.DefineBooleanMatrices(model); // optimization in this case
            else lagrangeEnumerator.DefineLagrangesAndBooleanMatrices(model);

            // Log dof statistics
            watch.Stop();
            Logger.LogTaskDuration("Dof ordering", watch.ElapsedMilliseconds);
            Logger.LogNumDofs("Global dofs", globalOrdering.NumGlobalFreeDofs);
            int numExpandedDomainFreeDofs = 0;
            foreach (var subdomain in model.Subdomains)
            {
                numExpandedDomainFreeDofs += subdomain.FreeDofOrdering.NumFreeDofs;
            }
            Logger.LogNumDofs("Expanded domain dofs", numExpandedDomainFreeDofs);
            Logger.LogNumDofs("Lagrange multipliers", lagrangeEnumerator.NumLagrangeMultipliers);

            //Leftover code from Model.ConnectDataStructures().
            //EnumerateSubdomainLagranges();
            //EnumerateDOFMultiplicity();
        }

        public void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        public void Solve()
        {
            var watch = new Stopwatch();
            foreach (ISingleSubdomainLinearSystem linearSystem in linearSystems.Values)
            {
                if (linearSystem.Solution == null) linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
            }

            // Calculate generalized inverses and rigid body modes of subdomains to assemble the interface flexibility matrix. 
            if (isStiffnessModified)
            {
                // Calculate the preconditioner before factorizing each subdomain's Kff 
                watch.Start();
                preconditioner = preconditionerFactory.CreatePreconditioner(stiffnessDistribution, dofSeparator,
                    lagrangeEnumerator, matrixManagersGeneral);
                watch.Stop();
                Logger.LogTaskDuration("Calculating preconditioner", watch.ElapsedMilliseconds);

                // Factorize each subdomain's Kff
                watch.Restart();
                foreach (int s in subdomains.Keys)
                {
                    matrixManagers[s].InvertKff(factorPivotTolerances[s], factorizeInPlace);
                }
                watch.Stop();
                Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);

                watch.Restart();
                BuildProjection();
                flexibility = new Feti1FlexibilityMatrix(matrixManagers, lagrangeEnumerator);
                watch.Stop();
                Logger.LogTaskDuration("Setting up interface problem", watch.ElapsedMilliseconds);
                watch.Reset();

                isStiffnessModified = false;
            }

            // Calculate the rhs vectors of the interface system
            watch.Start();
            Vector disconnectedDisplacements = CalcDisconnectedDisplacements();
            Vector rbmWork = CalcRigidBodyModesWork();
            double globalForcesNorm = CalcGlobalForcesNorm();
            watch.Stop();
            Logger.LogTaskDuration("Setting up interface problem", watch.ElapsedMilliseconds);

            // Solve the interface problem
            watch.Restart();
            Vector lagranges = interfaceProblemSolver.CalcLagrangeMultipliers(flexibility, preconditioner, projection,
                disconnectedDisplacements, rbmWork, globalForcesNorm, Logger);
            watch.Stop();
            Logger.LogTaskDuration("Solving interface problem", watch.ElapsedMilliseconds);

            // Calculate the displacements of each subdomain
            watch.Restart();
            Vector rbmCoeffs = CalcRigidBodyModesCoefficients(disconnectedDisplacements, lagranges);
            Dictionary<int, Vector> actualDisplacements = CalcActualDisplacements(lagranges, rbmCoeffs);
            foreach (var idSystem in linearSystems) idSystem.Value.SolutionConcrete = actualDisplacements[idSystem.Key];
            watch.Stop();
            Logger.LogTaskDuration("Calculate displacements from lagrange multipliers", watch.ElapsedMilliseconds);

            Logger.IncrementAnalysisStep();
        }

        /// <summary>
        /// Does not mutate this object.
        /// </summary>
        internal Dictionary<int, Vector> CalcActualDisplacements(Vector lagranges, Vector rigidBodyModeCoeffs)
        {
            var actualdisplacements = new Dictionary<int, Vector>();
            int rbmOffset = 0; //TODO: For this to work in parallel, each subdomain must store its offset.
            foreach (var linearSystem in linearSystems.Values)
            {
                // u = inv(K) * (f - B^T * λ), for non floating subdomains
                // u = generalizedInverse(K) * (f - B^T * λ) + R * a, for floating subdomains
                int s = linearSystem.Subdomain.ID;
                Vector forces = linearSystem.RhsConcrete - lagrangeEnumerator.BooleanMatrices[s].Multiply(lagranges, true);
                Vector displacements = matrixManagers[s].MultiplyInverseKffTimes(forces);

                foreach (Vector rbm in matrixManagers[s].RigidBodyModes)
                {
                    displacements.AxpyIntoThis(rbm, rigidBodyModeCoeffs[rbmOffset++]);
                }

                actualdisplacements[s] = displacements;
            }
            return actualdisplacements;
        }

        /// <summary>
        /// d = sum(Bs * generalInverse(Ks) * fs), where fs are the nodal forces applied to the dofs of subdomain s.
        /// Does not mutate this object.
        /// </summary>
        internal Vector CalcDisconnectedDisplacements()
        {
            var displacements = Vector.CreateZero(lagrangeEnumerator.NumLagrangeMultipliers);
            foreach (var linearSystem in linearSystems.Values)
            {
                int s = linearSystem.Subdomain.ID;
                Vector f = linearSystem.RhsConcrete;
                SignedBooleanMatrixColMajor B = lagrangeEnumerator.BooleanMatrices[s];
                Vector Kf = matrixManagers[s].MultiplyInverseKffTimes(f);
                Vector BKf = B.Multiply(Kf, false);
                displacements.AddIntoThis(BKf);
            }
            return displacements;
        }

        /// <summary>
        /// Does not mutate this object.
        /// </summary>
        internal Vector CalcRigidBodyModesCoefficients(Vector disconnectedDisplacements, Vector lagrangeMultipliers)
        {
            var flexibilityTimesLagranges = Vector.CreateZero(lagrangeEnumerator.NumLagrangeMultipliers);
            flexibility.Multiply(lagrangeMultipliers, flexibilityTimesLagranges);
            return projection.CalcRigidBodyModesCoefficients(flexibilityTimesLagranges, disconnectedDisplacements);
        }

        /// <summary>
        /// e = [ R1^T * f1; R2^T * f2; ... Rns^T * fns] 
        /// It doesn't mutate this object.
        /// </summary>
        internal Vector CalcRigidBodyModesWork() //TODO: this should probably be decomposed to subdomains
        {
            int workLength = 0;
            foreach (var linearSystem in linearSystems.Values)
            {
                workLength += matrixManagers[linearSystem.Subdomain.ID].RigidBodyModes.Count;
            }
            var work = Vector.CreateZero(workLength);

            int idx = 0;
            foreach (var linearSystem in linearSystems.Values)
            {
                int s = linearSystem.Subdomain.ID;
                Vector forces = linearSystem.RhsConcrete;
                foreach (Vector rbm in matrixManagers[s].RigidBodyModes) work[idx++] = rbm * forces;
            }

            return work;
        }

        private void BuildProjection()
        {
            if (!projectionMatrixQIsIdentity) // Previously (pde == PdeOrder.Fourth) || (!problemIsHomogeneous)
            {
                // Q = preconditioner
                projection = new Feti1Projection(lagrangeEnumerator.BooleanMatrices, matrixManagers,
                    new PreconditionerAsMatrixQ(preconditioner));
            }
            else
            {
                // Q = indentity
                projection = new Feti1Projection(lagrangeEnumerator.BooleanMatrices, matrixManagers, new IdentityMatrixQ());
            }
            projection.InvertCoarseProblemMatrix();
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
            subdomainGlobalMapping = new Feti1SubdomainGlobalMapping(model, dofSeparator, stiffnessDistribution);
        }

        public class Builder
        {
            private readonly Dictionary<int, double> factorPivotTolerances;
            private readonly IFeti1SubdomainMatrixManagerFactory matrixManagerFactory;

            public Builder(IFeti1SubdomainMatrixManagerFactory matrixManagerFactory, 
                Dictionary<int, double> factorPivotTolerances)
            {
                //TODO: This is a very volatile parameter and the user should not have to specify it. 
                this.factorPivotTolerances = factorPivotTolerances;
                this.matrixManagerFactory = matrixManagerFactory;
            }

            public IDofOrderer DofOrderer { get; set; } = 
                new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public IFeti1InterfaceProblemSolver InterfaceProblemSolver { get; set; } = 
                (new Feti1ProjectedInterfaceProblemSolver.Builder()).Build();

            public IFetiPreconditionerFactory PreconditionerFactory { get; set; } = new LumpedPreconditioner.Factory();
            public bool ProblemIsHomogeneous { get; set; } = true;
            public bool ProjectionMatrixQIsIdentity { get; set; } = true;
            //public PdeOrder PdeOrder { get; set; } = PdeOrder.Second; // Instead the user explicitly sets Q.

            public Feti1Solver BuildSolver(IStructuralModel model)
                => new Feti1Solver(model, matrixManagerFactory, DofOrderer, factorPivotTolerances, PreconditionerFactory,
                     InterfaceProblemSolver, ProblemIsHomogeneous, ProjectionMatrixQIsIdentity);
        }
    }
}
