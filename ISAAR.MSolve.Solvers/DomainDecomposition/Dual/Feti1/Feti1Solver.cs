using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Projection;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Rigid body modes do not have to be computed each time the stiffness matrix changes. 
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1
{
    public class Feti1Solver : ISolver_v2
    {
        internal const string name = "FETI Level-1 Solver"; // for error messages
        private readonly Dictionary<int, SkylineAssembler> assemblers;
        private readonly LagrangeMultipliersEnumerator lagrangeEnumerator;
        private readonly ICrosspointStrategy crosspointStrategy = new FullyRedundantConstraints();
        private readonly IDofOrderer dofOrderer;
        private readonly double factorizationPivotTolerance;
        private readonly IFeti1InterfaceProblemSolver interfaceProblemSolver;
        private readonly bool problemIsHomogeneous;
        private readonly bool projectionMatrixQIsIdentity;
        private readonly Dictionary<int, SingleSubdomainSystem<SkylineMatrix>> linearSystems;
        private readonly IStructuralModel_v2 model;
        //private readonly PdeOrder pde; // Instead the user explicitly sets Q.
        private readonly IFetiPreconditionerFactory preconditionerFactory;

        //TODO: fix the mess of Dictionary<int, ISubdomain>, List<ISubdomain>, Dictionary<int, Subdomain>, List<Subdomain>
        //      The concrete are useful for the preprocessor mostly, while analyzers, solvers need the interface versions.
        //      Lists are better in analyzers and solvers. Dictionaries may be better in the preprocessors.
        private readonly Dictionary<int, ISubdomain_v2> subdomains;

        private Feti1DofSeparator dofSeparator;
        private bool factorizeInPlace = true;
        private Dictionary<int, SemidefiniteCholeskySkyline> factorizations;
        private Feti1FlexibilityMatrix flexibility;
        private bool isStiffnessModified = true;
        private Dictionary<int, List<Vector>> rigidBodyModes;
        private IFetiPreconditioner preconditioner;
        private Feti1Projection projection;
        private IFeti1StiffnessDistribution stiffnessDistribution;

        private Feti1Solver(IStructuralModel_v2 model, IDofOrderer dofOrderer, double factorizationPivotTolerance,
            IFetiPreconditionerFactory preconditionerFactory, IFeti1InterfaceProblemSolver interfaceProblemSolver,
            bool problemIsHomogeneous, bool projectionMatrixQIsIdentity)
        {
            // Model
            if (model.Subdomains.Count == 1) throw new InvalidSolverException(
                $"{name} cannot be used if there is only 1 subdomain");
            this.model = model;

            // Subdomains
            subdomains = new Dictionary<int, ISubdomain_v2>();
            foreach (Subdomain_v2 subdomain in model.Subdomains) subdomains[subdomain.ID] = subdomain;

            // Linear systems
            linearSystems = new Dictionary<int, SingleSubdomainSystem<SkylineMatrix>>();
            var tempLinearSystems = new Dictionary<int, ILinearSystem_v2>();
            assemblers = new Dictionary<int, SkylineAssembler>();
            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                var linearSystem = new SingleSubdomainSystem<SkylineMatrix>(subdomain);
                linearSystems.Add(id, linearSystem);
                tempLinearSystems.Add(id, linearSystem);
                linearSystem.MatrixObservers.Add(this);
                assemblers.Add(id, new SkylineAssembler());
            }
            LinearSystems = tempLinearSystems;

            this.dofOrderer = dofOrderer;
            this.factorizationPivotTolerance = factorizationPivotTolerance;
            this.preconditionerFactory = preconditionerFactory;

            // PCPG
            this.interfaceProblemSolver = interfaceProblemSolver;

            this.lagrangeEnumerator = new LagrangeMultipliersEnumerator(crosspointStrategy);

            // Homogeneous/heterogeneous problems
            this.problemIsHomogeneous = problemIsHomogeneous;
            this.projectionMatrixQIsIdentity = projectionMatrixQIsIdentity;
        }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }
        public DualSolverLogger Logger { get; } = new DualSolverLogger();

        public Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider_v2 elementMatrixProvider)
        {
            var matricesInternal = new Dictionary<int, IMatrixView>();
            var matricesResult = new Dictionary<int, IMatrix>();
            foreach (ISubdomain_v2 subdomain in model.Subdomains) //TODO: this must be done in parallel
            {
                SkylineMatrix matrix = assemblers[subdomain.ID].BuildGlobalMatrix(subdomain.FreeDofOrdering,
                    subdomain.Elements, elementMatrixProvider);
                matricesInternal[subdomain.ID] = matrix;
                matricesResult[subdomain.ID] = matrix;
            }

            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            DetermineStiffnessDistribution(matricesInternal);

            return matricesResult;
        }

        public Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr)> BuildGlobalSubmatrices(IElementMatrixProvider_v2 elementMatrixProvider)
        {
            var matricesResult = new Dictionary<int, (IMatrix Aff, IMatrixView Afc, IMatrixView Acf, IMatrixView Acc)>();
            var matricesInternal = new Dictionary<int, IMatrixView>();
            foreach (ISubdomain_v2 subdomain in model.Subdomains) //TODO: this must be done in parallel
            {
                if (subdomain.ConstrainedDofOrdering == null)
                {
                    throw new InvalidOperationException("In order to build the matrices corresponding to constrained dofs of,"
                        + $" subdomain {subdomain.ID}, they must have been ordered first.");
                }
                (SkylineMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) =
                    assemblers[subdomain.ID].BuildGlobalSubmatrices(subdomain.FreeDofOrdering,
                    subdomain.ConstrainedDofOrdering, subdomain.Elements, elementMatrixProvider);
                matricesResult[subdomain.ID] = (Kff, Kfc, Kcf, Kcc);
                matricesInternal[subdomain.ID] = Kff;
            }

            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            DetermineStiffnessDistribution(matricesInternal);

            return matricesResult;
        }

        //TODO: this and the fields should be handled by a class that handles dof mappings.
        public Dictionary<int, SparseVector> DistributeNodalLoads(Table<INode, DOFType, double> globalNodalLoads)
            => stiffnessDistribution.SubdomainGlobalConversion.DistributeNodalLoads(LinearSystems, globalNodalLoads);

        //TODO: this and the fields should be handled by a class that handles dof mappings.
        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainDisplacements)
            => stiffnessDistribution.SubdomainGlobalConversion.GatherGlobalDisplacements(subdomainDisplacements);

        public void HandleMatrixWillBeSet()
        {
            isStiffnessModified = true;
            factorizations = null;
            flexibility = null;
            rigidBodyModes = null;
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
            // Order dofs
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderFreeDofs(model);
            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                assemblers[subdomain.ID].HandleDofOrderingWillBeModified();
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
                if (alsoOrderConstrainedDofs) subdomain.ConstrainedDofOrdering = dofOrderer.OrderConstrainedDofs(subdomain);

                // The next must done by the analyzer, so that subdomain.Forces is retained when doing back to back analyses.
                //subdomain.Forces = linearSystem.CreateZeroVector();
            }

            // Define boundary / internal dofs
            dofSeparator = new Feti1DofSeparator();
            dofSeparator.SeparateBoundaryInternalDofs(model);

            // Define lagrange multipliers and boolean matrices
            if (problemIsHomogeneous) lagrangeEnumerator.DefineBooleanMatrices(model, dofSeparator); // optimization in this case
            else lagrangeEnumerator.DefineLagrangesAndBooleanMatrices(model, dofSeparator);

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
            }

            //Leftover code from Model.ConnectDataStructures().
            //EnumerateSubdomainLagranges();
            //EnumerateDOFMultiplicity();
        }

        public void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        public void Solve()
        {
            foreach (var linearSystem in linearSystems.Values)
            {
                if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            }

            // Calculate generalized inverses and rigid body modes of subdomains to assemble the interface flexibility matrix. 
            if (isStiffnessModified)
            {
                BuildPreconditioner();
                FactorizeMatrices();
                BuildProjection();
                isStiffnessModified = false;
                flexibility = new Feti1FlexibilityMatrix(factorizations, lagrangeEnumerator);
            }

            // Calculate the rhs vectors of the interface system
            Vector disconnectedDisplacements = CalcDisconnectedDisplacements();
            Vector rbmWork = CalcRigidBodyModesWork();
            double globalForcesNorm = CalcGlobalForcesNorm();

            // Solve the interface problem
            Vector lagranges = interfaceProblemSolver.CalcLagrangeMultipliers(flexibility, preconditioner, projection,
                disconnectedDisplacements, rbmWork, globalForcesNorm, Logger);

            // Calculate the displacements of each subdomain
            Vector rbmCoeffs = CalcRigidBodyModesCoefficients(disconnectedDisplacements, lagranges);
            Dictionary<int, Vector> actualDisplacements = CalcActualDisplacements(lagranges, rbmCoeffs);
            foreach (var idSystem in linearSystems) idSystem.Value.Solution = actualDisplacements[idSystem.Key];
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
                int id = linearSystem.Subdomain.ID;
                Vector forces = linearSystem.RhsVector - lagrangeEnumerator.BooleanMatrices[id].Multiply(lagranges, true);
                Vector displacements = factorizations[id].MultiplyGeneralizedInverseMatrixTimesVector(forces);

                foreach (Vector rbm in rigidBodyModes[id]) displacements.AxpyIntoThis(rbm, rigidBodyModeCoeffs[rbmOffset++]);

                actualdisplacements[id] = displacements;
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
                int id = linearSystem.Subdomain.ID;
                Vector f = linearSystem.RhsVector;
                SignedBooleanMatrix boolean = lagrangeEnumerator.BooleanMatrices[id];
                Vector Kf = factorizations[id].MultiplyGeneralizedInverseMatrixTimesVector(f);
                Vector BKf = boolean.Multiply(Kf, false);
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
                workLength += rigidBodyModes[linearSystem.Subdomain.ID].Count;
            }
            var work = Vector.CreateZero(workLength);

            int idx = 0;
            foreach (var linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                Vector forces = linearSystem.RhsVector;
                foreach (Vector rbm in rigidBodyModes[id]) work[idx++] = rbm * forces;
            }

            return work;
        }

        private void BuildPreconditioner()
        {
            // Create the preconditioner. 
            //TODO: this should be done simultaneously with the factorizations to avoid duplicate factorizations.
            var stiffnessMatrices = new Dictionary<int, IMatrixView>();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                stiffnessMatrices.Add(linearSystem.Subdomain.ID, linearSystem.Matrix);
            }
            this.preconditioner = preconditionerFactory.CreatePreconditioner(stiffnessDistribution, dofSeparator,
                lagrangeEnumerator, stiffnessMatrices);
        }

        private void BuildProjection()
        {
            if (!projectionMatrixQIsIdentity) // Previously (pde == PdeOrder.Fourth) || (!problemIsHomogeneous)
            {
                // Q = preconditioner
                projection = new Feti1Projection(lagrangeEnumerator.BooleanMatrices, rigidBodyModes,
                    new PreconditionerAsMatrixQ(preconditioner));
            }
            else
            {
                // Q = indentity
                projection = new Feti1Projection(lagrangeEnumerator.BooleanMatrices, rigidBodyModes, new IdentityMatrixQ());
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
                subdomainForces[linearSystem.Subdomain.ID] = linearSystem.RhsVector;
            }
            return stiffnessDistribution.SubdomainGlobalConversion.CalculateGlobalForcesNorm(subdomainForces);
        }

        private void DetermineStiffnessDistribution(Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            if (problemIsHomogeneous) stiffnessDistribution = new HomogeneousStiffnessDistribution(model, dofSeparator);
            else stiffnessDistribution = new HeterogeneousStiffnessDistribution(model, dofSeparator, stiffnessMatrices);
        }

        private void FactorizeMatrices()
        {
            factorizations = new Dictionary<int, SemidefiniteCholeskySkyline>();
            rigidBodyModes = new Dictionary<int, List<Vector>>();
            foreach (var linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                factorizations[id] =
                    linearSystem.Matrix.FactorSemidefiniteCholesky(factorizeInPlace, factorizationPivotTolerance);
                rigidBodyModes[id] = new List<Vector>();
                foreach (double[] rbm in factorizations[id].NullSpaceBasis)
                {
                    rigidBodyModes[id].Add(Vector.CreateFromArray(rbm, false));
                }
            }
        }

        public class Builder
        {
            private readonly double factorizationPivotTolerance;

            public Builder(double factorizationPivotTolerance)
            {
                //TODO: This is a very volatile parameter and the user should not have to specify it. 
                this.factorizationPivotTolerance = factorizationPivotTolerance;
            }

            public IDofOrderer DofOrderer { get; set; } = 
                new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public IFeti1InterfaceProblemSolver InterfaceProblemSolver { get; set; } = 
                (new Feti1ProjectedInterfaceProblemSolver.Builder()).Build();

            public IFetiPreconditionerFactory PreconditionerFactory { get; set; } = new Feti1LumpedPreconditioner.Factory();
            public bool ProblemIsHomogeneous { get; set; } = true;
            public bool ProjectionMatrixQIsIdentity { get; set; } = true;
            //public PdeOrder PdeOrder { get; set; } = PdeOrder.Second; // Instead the user explicitly sets Q.

            public Feti1Solver BuildSolver(Model_v2 model)
                => new Feti1Solver(model, DofOrderer, factorizationPivotTolerance, PreconditionerFactory,
                     InterfaceProblemSolver, ProblemIsHomogeneous, ProjectionMatrixQIsIdentity);
        }
    }
}
