using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Rigid body modes do not have to be computed each time the stiffness matrix changes. 
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Feti1
{
    public class Feti1Solver : ISolver_v2
    {
        internal const string name = "FETI-1 Solver"; // for error messages
        private readonly Dictionary<int, SkylineAssembler> assemblers;
        private readonly LagrangeMultipliersEnumerator lagrangeEnumerator;
        private readonly ICrosspointStrategy crosspointStrategy = new FullyRedundantConstraints();
        private readonly IDofOrderer dofOrderer;
        private readonly double factorizationPivotTolerance;
        private readonly IFeti1InterfaceProblemSolver interfaceProblemSolver;
        private readonly bool isProblemHomogeneous;
        private readonly Dictionary<int, SingleSubdomainSystem<SkylineMatrix>> linearSystems;
        private readonly IStructuralModel_v2 model;
        private readonly IExactResidualCalculator pcpgExactResidual; //TODO: this and related logic should be move to the class that handles PCG/PCPG
        private readonly PdeOrder pde;
        private readonly IFetiPreconditionerFactory preconditionerFactory;

        //TODO: fix the mess of Dictionary<int, ISubdomain>, List<ISubdomain>, Dictionary<int, Subdomain>, List<Subdomain>
        //      The concrete are useful for the preprocessor mostly, while analyzers, solvers need the interface versions.
        //      Lists are better in analyzers and solvers. Dictionaries may be better in the preprocessors.
        private readonly Dictionary<int, ISubdomain_v2> subdomains;

        private DofSeparator dofSeparator;
        private bool factorizeInPlace = true;
        private Dictionary<int, SemidefiniteCholeskySkyline> factorizations;
        private bool isStiffnessModified = true;
        private Dictionary<int, List<Vector>> rigidBodyModes;
        private IFetiPreconditioner preconditioner;
        private Feti1Projection projection;
        private IStiffnessDistribution stiffnessDistribution;

        private Feti1Solver(IStructuralModel_v2 model, IDofOrderer dofOrderer, double factorizationPivotTolerance,
            IFetiPreconditionerFactory preconditionerFactory, IFeti1InterfaceProblemSolver interfaceProblemSolver,
            IExactResidualCalculator pcgExactResidual, bool isProblemHomogeneous, PdeOrder pde)
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
            this.pcpgExactResidual = pcgExactResidual;

            this.lagrangeEnumerator = new LagrangeMultipliersEnumerator(crosspointStrategy);

            // Homogeneous/heterogeneous problems
            this.pde = pde;
            this.isProblemHomogeneous = isProblemHomogeneous;
        }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }
        public FetiLogger Logger { get; } = new FetiLogger();

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
            dofSeparator = new DofSeparator();
            dofSeparator.SeparateBoundaryInternalDofs(model);

            // Define lagrange multipliers and boolean matrices
            if (isProblemHomogeneous) lagrangeEnumerator.DefineBooleanMatrices(model, dofSeparator); // optimization in this case
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
                preconditioner = CalcPreconditioner();
                FactorizeMatrices();
                projection = CalcProjection(preconditioner);
                isStiffnessModified = false;
            }
            var flexibility = new Feti1FlexibilityMatrix(factorizations, lagrangeEnumerator);

            // Calculate the rhs vectors of the interface system
            Vector disconnectedDisplacements = CalcDisconnectedDisplacements();
            Vector rbmWork = CalcRigidBodyModesWork();
            double globalForcesNorm = CalcGlobalForcesNorm();

            //TODO: This should be done in the method that handles PCG/PCPG
            //// Handle exact residual calculations if they are enabled during testing
            //CalculateExactResidualNorm calcExactResidualNorm = null;
            //if (pcpgExactResidual != null)
            //{
            //    calcExactResidualNorm = lagr => CalculateExactResidualNorm(lagr, flexibility, projection,
            //        disconnectedDisplacements);
            //}

            // Solve the interface problem
            Vector lagranges = interfaceProblemSolver.CalcLagrangeMultipliers(flexibility, preconditioner, projection,
                disconnectedDisplacements, rbmWork, globalForcesNorm, Logger);

            // Calculate the displacements of each subdomain
            Vector rbmCoeffs = CalcRigidBodyModesCoefficients(flexibility, projection, disconnectedDisplacements,
                lagranges);
            Dictionary<int, Vector> actualDisplacements = CalcActualDisplacements(lagranges, rbmCoeffs);
            foreach (var idSystem in linearSystems) idSystem.Value.Solution = actualDisplacements[idSystem.Key];
        }

        private Dictionary<int, Vector> CalcActualDisplacements(Vector lagranges, Vector rigidBodyModeCoeffs)
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
        /// </summary>
        private Vector CalcDisconnectedDisplacements()
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

        //TODO: This is only used by the corresponding PCPG convergence criteria, which itself is for testing purposes only. 
        //      It should be moved to another class. If private members of this class are needed, then reflection should be used
        //      to access them.
        private double CalcExactResidualNorm(Vector lagranges, Feti1FlexibilityMatrix flexibility, 
            Feti1Projection projection, Vector disconnectedDisplacements)
        {
            Vector rbmCoeffs = CalcRigidBodyModesCoefficients(flexibility, projection, disconnectedDisplacements, lagranges);
            var subdomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var idDisplacements in CalcActualDisplacements(lagranges, rbmCoeffs))
            {
                subdomainDisplacements[idDisplacements.Key] = idDisplacements.Value;
            }
            Vector globalDisplacements = GatherGlobalDisplacements(subdomainDisplacements);
            return pcpgExactResidual.CalculateExactResidualNorm(globalDisplacements);
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

        private IFetiPreconditioner CalcPreconditioner()
        {
            // Create the preconditioner. 
            //TODO: this should be done simultaneously with the factorizations to avoid duplicate factorizations.
            var stiffnessMatrices = new Dictionary<int, IMatrixView>();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                stiffnessMatrices.Add(linearSystem.Subdomain.ID, linearSystem.Matrix);
            }
            return preconditionerFactory.CreatePreconditioner(stiffnessDistribution, dofSeparator, lagrangeEnumerator,
                stiffnessMatrices);
        }

        private Feti1Projection CalcProjection(IFetiPreconditioner preconditioner)
        {
            Feti1Projection projection;
            if ((pde == PdeOrder.Fourth) || (!isProblemHomogeneous))
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
            return projection;
        }

        private Vector CalcRigidBodyModesCoefficients(Feti1FlexibilityMatrix flexibility, Feti1Projection projection,
            Vector disconnectedDisplacements, Vector lagrangeMultipliers)
        {
            var flexibilityTimesLagranges = Vector.CreateZero(lagrangeEnumerator.NumLagrangeMultipliers);
            flexibility.Multiply(lagrangeMultipliers, flexibilityTimesLagranges);
            return projection.CalcRigidBodyModesCoefficients(flexibilityTimesLagranges, disconnectedDisplacements);
        }

        /// <summary>
        /// e = [ R1^T * f1; R2^T * f2; ... Rns^T * fns] 
        /// </summary>
        private Vector CalcRigidBodyModesWork() //TODO: this should probably be decomposed to subdomains
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

        private void DetermineStiffnessDistribution(Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            if (isProblemHomogeneous) stiffnessDistribution = new HomogeneousStiffnessDistribution(model, dofSeparator);
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
                new Feti1ProjectedInterfaceProblemSolver(1E-7, 1.0, Feti1ProjectedInterfaceProblemSolver.ProjectionSide.Both,
                    Feti1ProjectedInterfaceProblemSolver.ProjectionSide.Both);

            public bool IsProblemHomogeneous { get; set; } = true;
            internal IExactResidualCalculator PcgExactResidual { get; set; } = null;
            public double PcgConvergenceTolerance { get; set; } = 1E-7;
            public double PcgMaxIterationsOverSize { get; set; } = 1.0;
            public PdeOrder PdeOrder { get; set; } = PdeOrder.Second;
            public IFetiPreconditionerFactory PreconditionerFactory { get; set; } = new Feti1LumpedPreconditioner.Factory();
            
            public Feti1Solver BuildSolver(Model_v2 model)
                => new Feti1Solver(model, DofOrderer, factorizationPivotTolerance, PreconditionerFactory,
                     InterfaceProblemSolver, PcgExactResidual, IsProblemHomogeneous, PdeOrder);
        }
    }
}
