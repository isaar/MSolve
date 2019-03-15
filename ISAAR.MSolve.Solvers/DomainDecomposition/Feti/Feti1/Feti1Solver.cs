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
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Rigid body modes do not have to be computed each time the stiffness matrix changes. 
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Feti1
{
    public class Feti1Solver : ISolver_v2
    {
        private readonly Dictionary<int, SkylineAssembler> assemblers;
        private readonly LagrangeMultipliersEnumerator lagrangeEnumerator;
        private readonly ICrosspointStrategy crosspointStrategy = new FullyRedundantConstraints();
        private readonly IDofOrderer dofOrderer;
        private readonly double factorizationPivotTolerance;
        private readonly bool isProblemHomogeneous;
        private readonly Dictionary<int, SingleSubdomainSystem<SkylineMatrix>> linearSystems;
        private readonly FetiLogger logger;
        private readonly IStructuralModel_v2 model;
        private readonly string name = "FETI-1 Solver"; // for error messages
        private readonly IExactResidualCalculator pcpgExactResidual;
        private readonly double pcpgConvergenceTolerance;
        private readonly double pcpgMaxIterationsOverSize;
        private readonly IFetiPreconditionerFactory preconditionerFactory;
        private readonly IStiffnessDistribution stiffnessDistribution;

        //TODO: fix the mess of Dictionary<int, ISubdomain>, List<ISubdomain>, Dictionary<int, Subdomain>, List<Subdomain>
        //      The concrete are useful for the preprocessor mostly, while analyzers, solvers need the interface versions.
        //      Lists are better in analyzers and solvers. Dictionaries may be better in the preprocessors.
        private readonly Dictionary<int, ISubdomain_v2> subdomains; 

        private DofSeparator dofSeparator;
        private bool factorizeInPlace = true;
        private bool mustFactorize = true;
        private Dictionary<int, SemidefiniteCholeskySkyline> factorizations;
        private Dictionary<int, List<Vector>> rigidBodyModes;

        private Feti1Solver(IStructuralModel_v2 model, IDofOrderer dofOrderer, double factorizationPivotTolerance, 
            double pcpgMaxIterationsOverSize, double pcpgConvergenceTolerance, IExactResidualCalculator pcpgExactResidual,
            IFetiPreconditionerFactory preconditionerFactory, bool isProblemHomogeneous, FetiLogger logger)
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
            this.logger = logger;

            // PCPG
            this.pcpgMaxIterationsOverSize = pcpgMaxIterationsOverSize;
            this.pcpgConvergenceTolerance = pcpgConvergenceTolerance;
            this.pcpgExactResidual = pcpgExactResidual;

            this.lagrangeEnumerator = new LagrangeMultipliersEnumerator(crosspointStrategy);

            // Homegeneous/heterogeneous problems
            this.isProblemHomogeneous = isProblemHomogeneous;
            if (isProblemHomogeneous) this.stiffnessDistribution = new HomogeneousStiffnessDistribution();
            else this.stiffnessDistribution = new HeterogeneousStiffnessDistribution();
        }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider_v2 elementMatrixProvider)
        {
            return assemblers[subdomain.ID].BuildGlobalMatrix(
                subdomain.FreeDofOrdering, subdomain.Elements, elementMatrixProvider);
        }

        public (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree, 
            IMatrixView matrixConstrConstr) BuildGlobalSubmatrices(
            ISubdomain_v2 subdomain, IElementMatrixProvider_v2 elementMatrixProvider)
        {
            if (subdomain.ConstrainedDofOrdering == null)
            {
                throw new InvalidOperationException("In order to build the matrices corresponding to constrained dofs,"
                    + " they must have been ordered first.");
            }
            return assemblers[subdomain.ID].BuildGlobalSubmatrices(subdomain.FreeDofOrdering, subdomain.ConstrainedDofOrdering,
                subdomain.Elements, elementMatrixProvider);
        }

        //TODO: this is for homogeneous problems only.
        public Dictionary<int, SparseVector> DistributeNodalLoads(Table<INode, DOFType, double> globalNodalLoads)
            => stiffnessDistribution.NodalLoadDistributor.DistributeNodalLoads(LinearSystems, globalNodalLoads);

        //TODO: this and the fields should be handled by a class that handles dofs.
        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainDisplacements) 
        {
            var globalDisplacements = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            foreach (var subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                IVectorView displacements = subdomainDisplacements[id]; //TODO: benchmark the performance if this was concrete Vector

                // Internal dofs are copied as is.
                foreach (int internalDof in dofSeparator.InternalDofs[id])
                {
                    int globalDofIdx = subdomainToGlobalDofs[internalDof];
                    globalDisplacements[globalDofIdx] = displacements[internalDof];
                }

                // For boundary dofs we take the mean value across subdomains. 
                for (int i = 0; i < dofSeparator.BoundaryDofs[id].Length; ++i)
                {
                    int multiplicity = dofSeparator.BoundaryDofsMultiplicity[id][i];
                    int subdomainDofIdx = dofSeparator.BoundaryDofs[id][i];
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                    globalDisplacements[globalDofIdx] += displacements[subdomainDofIdx] / multiplicity;
                }
            }
            return globalDisplacements;
        }

        public void HandleMatrixWillBeSet()
        {
            mustFactorize = true;
            factorizations = null;
        }

        public void Initialize()
        {
            
        }

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
            if (isProblemHomogeneous) lagrangeEnumerator.DefineBooleanMatrices(model, dofSeparator);
            else lagrangeEnumerator.DefineLagrangesAndBooleanMatrices(model, dofSeparator);

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

            // Create the preconditioner. 
            //TODO: this should be done simultaneously with the factorizations to avoid duplicate factorizations.
            var stiffnessMatrices = new Dictionary<int, IMatrixView>();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                stiffnessMatrices.Add(linearSystem.Subdomain.ID, linearSystem.Matrix);
            }
            if (!isProblemHomogeneous)
            {
                Table<INode, DOFType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses
                    = BoundaryDofLumpedStiffness.ExtractBoundaryDofLumpedStiffnesses(dofSeparator, stiffnessMatrices);
                stiffnessDistribution.StoreStiffnesses(stiffnessMatrices, boundaryDofStiffnesses);
            }
            IFetiPreconditioner preconditioner = preconditionerFactory.CreatePreconditioner(stiffnessDistribution, dofSeparator, 
                lagrangeEnumerator, stiffnessMatrices);

            // Calculate generalized inverses and rigid body modes of subdomains to assemble the interface flexibility matrix. 
            if (mustFactorize)
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
                mustFactorize = false;
            }
            var flexibility = new Feti1FlexibilityMatrix(factorizations, lagrangeEnumerator);

            // Calculate the rhs vectors of the interface system
            Vector disconnectedDisplacements = CalculateDisconnectedDisplacements();
            Vector rbmWork = CalculateRigidBodyModesWork();

            // Define and initilize the projection
            var Q = Matrix.CreateIdentity(lagrangeEnumerator.NumLagrangeMultipliers);
            var projection = new Feti1Projection(lagrangeEnumerator.BooleanMatrices, rigidBodyModes, Q);
            projection.InvertCoarseProblemMatrix();

            // Calculate the norm of the forces vector Ku=f
            //TODO: It would be better to do that using the global vector to avoid the homogeneous/heterogeneous averaging
            //      That would require the analyzer to build the global vector too though. Caution: we cannot take just 
            //      the nodal loads from the model, since the effect of other loads is only taken into account int 
            //      linearSystem.Rhs. Even if we could easily create the global forces vector, it might be wrong since 
            //      the analyzer may modify some of these loads, depending on time step, loading step, etc.
            var subdomainForces = new Dictionary<int, Vector>();
            foreach (var linearSystem in linearSystems.Values)
            {
                subdomainForces[linearSystem.Subdomain.ID] = linearSystem.RhsVector;
            }
            double globalForcesNorm = CalculateGlobalForcesNorm(subdomainForces);

            // Handle exact residual calculations if they are enabled during testing
            CalculateExactResidualNorm calcExactResidualNorm = null;
            if (pcpgExactResidual != null)
            {
                calcExactResidualNorm = lagranges => CalculateExactResidualNorm(lagranges, flexibility, projection, 
                    disconnectedDisplacements);
            }

            // Run the PCPG algorithm
            var pcpg = new PcpgAlgorithm(pcpgMaxIterationsOverSize, pcpgConvergenceTolerance, calcExactResidualNorm);
            var lagrangeMultipliers = Vector.CreateZero(lagrangeEnumerator.NumLagrangeMultipliers);
            PcpgStatistics stats = pcpg.Solve(flexibility, preconditioner, projection, disconnectedDisplacements, rbmWork,
                globalForcesNorm, lagrangeMultipliers);
            if (logger != null)
            {
                logger.NumUniqueGlobalFreeDofs = model.GlobalDofOrdering.NumGlobalFreeDofs;
                logger.NumExpandedDomainFreeDofs = 0;
                foreach (var subdomain in model.Subdomains)
                {
                    logger.NumExpandedDomainFreeDofs += subdomain.FreeDofOrdering.NumFreeDofs;
                }
                logger.NumLagrangeMultipliers = lagrangeMultipliers.Length;
                logger.PcpgIterations = stats.NumIterations;
                logger.PcpgResidualNormEstimateRatio = stats.ResidualNormEstimateRatio;
            }

            // Calculate the displacements of each subdomain
            Vector rbmCoeffs = CalculateRigidBodyModesCoefficients(flexibility, projection, disconnectedDisplacements, 
                lagrangeMultipliers);
            Dictionary<int, Vector> actualDisplacements = CalculateActualDisplacements(lagrangeMultipliers, rbmCoeffs);
            foreach (var idSystem in linearSystems) idSystem.Value.Solution = actualDisplacements[idSystem.Key];
        }

        private Dictionary<int, Vector> CalculateActualDisplacements(Vector lagranges, Vector rigidBodyModeCoeffs)
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
        private Vector CalculateDisconnectedDisplacements()
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
        /// This is only used by the corresponeing PCPG convergence criteria, which itself is fro testing purposes only.
        /// </summary>
        private double CalculateExactResidualNorm(Vector lagranges, Feti1FlexibilityMatrix flexibility, //TODO: move to another class
            Feti1Projection projection, Vector disconnectedDisplacements)
        {
            Vector rbmCoeffs = CalculateRigidBodyModesCoefficients(flexibility, projection, disconnectedDisplacements, lagranges);
            var subdomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var idDisplacements in CalculateActualDisplacements(lagranges, rbmCoeffs))
            {
                subdomainDisplacements[idDisplacements.Key] = idDisplacements.Value;
            }
            Vector globalDisplacements = GatherGlobalDisplacements(subdomainDisplacements);
            return pcpgExactResidual.CalculateExactResidualNorm(globalDisplacements);
        }

        //TODO: this should be used for non linear analyzers as well (instead of building the global RHS)
        //TODO: this should be implemented by a different class.
        //TODO: Is this correct? For the residual, it would be wrong to find f-K*u for each subdomain and then call this.
        private double CalculateGlobalForcesNorm(Dictionary<int, Vector> subdomainForces) //TODO: This should take into account the heterogenity.
        {
            double globalSum = 0.0;
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                double subdomainSum = 0.0;
                Vector forces = subdomainForces[id];
                foreach (int internalDof in dofSeparator.InternalDofs[id]) subdomainSum += forces[internalDof];
                for (int i = 0; i < dofSeparator.BoundaryDofs[id].Length; ++i)
                {
                    // E.g. 2 subdomains. Originally: sum += f * f. Now: sum += (2*f/2)(2*f/2)/2 + (2*f/2)(2*f/2)/2
                    // WARNING: This works only if nodal loads are distributed evenly among subdomains.
                    int multiplicity = dofSeparator.BoundaryDofsMultiplicity[id][i];
                    double totalForce = forces[dofSeparator.BoundaryDofs[id][i]] * multiplicity;
                    subdomainSum += (totalForce * totalForce) / multiplicity;
                }
                globalSum += subdomainSum;
            }
            return Math.Sqrt(globalSum);
        }

        private Vector CalculateRigidBodyModesCoefficients(Feti1FlexibilityMatrix flexibility, Feti1Projection projection,
            Vector disconnectedDisplacements, Vector lagrangeMultipliers)
        {
            var flexibilityTimesLagranges = Vector.CreateZero(lagrangeEnumerator.NumLagrangeMultipliers);
            flexibility.Multiply(lagrangeMultipliers, flexibilityTimesLagranges);
            return projection.CalculateRigidBodyModesCoefficients(flexibilityTimesLagranges, disconnectedDisplacements);
        }

        /// <summary>
        /// e = [ R1^T * f1; R2^T * f2; ... Rns^T * fns] 
        /// </summary>
        private Vector CalculateRigidBodyModesWork() //TODO: this should probably be decomposed to subdomains
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

        public class Builder
        {
            private readonly double factorizationPivotTolerance;

            public Builder(double factorizationPivotTolerance)
            {
                //TODO: This is a very volatile parameter and the user should not have to specify it. 
                this.factorizationPivotTolerance = factorizationPivotTolerance;
            }

            public IDofOrderer DofOrderer { get; set; } = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            public bool IsProblemHomogeneous { get; set; } = true;
            public FetiLogger Logger { get; set; } = null;
            internal IExactResidualCalculator PcpgExactResidual { get; set; } = null;
            public double PcpgConvergenceTolerance { get; set; } = 1E-7;
            public double PcpgMaxIterationsOverSize { get; set; } = 1.0;
            public IFetiPreconditionerFactory PreconditionerFactory { get; set; } = new Feti1LumpedPreconditioner.Factory();
            //public IStiffnessDistribution StiffnessDistribution { get; set; } = new HomogeneousStiffnessDistribution();

            public Feti1Solver BuildSolver(Model_v2 model)
                => new Feti1Solver(model, DofOrderer, factorizationPivotTolerance, PcpgMaxIterationsOverSize, 
                    PcpgConvergenceTolerance, PcpgExactResidual, PreconditionerFactory, IsProblemHomogeneous, Logger);
        }
    }
}
