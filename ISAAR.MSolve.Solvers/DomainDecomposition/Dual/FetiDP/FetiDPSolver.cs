using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
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
//TODO: Optimizations for the case that stiffness changes, but connectivity remains the same!
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    public class FetiDPSolver : ISolver
    {
        internal const string name = "FETI-DP Solver"; // for error messages
        private readonly IFetiDPCoarseProblemSolver coarseProblemSolver;
        private readonly ICornerNodeSelection cornerNodeSelection;
        private readonly ICrosspointStrategy crosspointStrategy = new FullyRedundantConstraints();
        private readonly IDofOrderer dofOrderer;
        private readonly FetiDPDofSeparator dofSeparator;
        private readonly IFetiDPInterfaceProblemSolver interfaceProblemSolver;
        private readonly Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers;
        private readonly Dictionary<int, IFetiSubdomainMatrixManager> matrixManagersGeneral; //TODO: redesign. They are the same as above, but Dictionary is not covariant
        private readonly Dictionary<int, ISingleSubdomainLinearSystem> linearSystems;
        private readonly IStructuralModel model;
        private readonly IFetiPreconditionerFactory preconditionerFactory;
        private readonly bool problemIsHomogeneous;
        private readonly IStiffnessDistribution stiffnessDistribution;

        //TODO: fix the mess of Dictionary<int, ISubdomain>, List<ISubdomain>, Dictionary<int, Subdomain>, List<Subdomain>
        //      The concrete are useful for the preprocessor mostly, while analyzers, solvers need the interface versions.
        //      Lists are better in analyzers and solvers. Dictionaries may be better in the preprocessors.
        private readonly Dictionary<int, ISubdomain> subdomains;

        private bool factorizeInPlace = true;
        private FetiDPFlexibilityMatrix flexibility;
        private bool isStiffnessModified = true;
        private FetiDPLagrangeMultipliersEnumerator lagrangeEnumerator;
        private IFetiPreconditioner preconditioner;
        private FetiDPSubdomainGlobalMapping subdomainGlobalMapping;

        private FetiDPSolver(IStructuralModel model, ICornerNodeSelection cornerNodeSelection,
            IFetiDPSubdomainMatrixManagerFactory matrixManagerFactory, IDofOrderer dofOrderer, 
            IFetiPreconditionerFactory preconditionerFactory, bool problemIsHomogeneous, 
            IFetiDPInterfaceProblemSolver interfaceProblemSolver)
        {
            // Model
            if (model.Subdomains.Count == 1) throw new InvalidSolverException(
                $"{name} cannot be used if there is only 1 subdomain");
            this.model = model;
            this.cornerNodeSelection = cornerNodeSelection;

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
                int s = subdomain.ID;
                var matrixManager = matrixManagerFactory.CreateMatricesManager(subdomain);
                matrixManagers[s] = matrixManager;
                matrixManagersGeneral[s] = matrixManager;
                this.linearSystems[s] = matrixManager.LinearSystem;
                externalLinearSystems[s] = matrixManager.LinearSystem;

                //TODO: This will call HandleMatrixWillBeSet() once for each subdomain. For now I will clear the data when 
                //      BuildMatrices() is called. Redesign this.
                //matrixManager.LinearSystem.MatrixObservers.Add(this); 
            }
            LinearSystems = externalLinearSystems;

            this.dofOrderer = dofOrderer;
            this.dofSeparator = new FetiDPDofSeparator();
            this.preconditionerFactory = preconditionerFactory;

            // Interface problem
            this.coarseProblemSolver = matrixManagerFactory.CreateCoarseProblemSolver(model.Subdomains);
            this.interfaceProblemSolver = interfaceProblemSolver;

            // Homogeneous/heterogeneous problems
            this.problemIsHomogeneous = problemIsHomogeneous;
            if (problemIsHomogeneous) this.stiffnessDistribution = new HomogeneousStiffnessDistribution(model, dofSeparator);
            else this.stiffnessDistribution = new HeterogeneousStiffnessDistribution(model, dofSeparator);
        }

        public Dictionary<int, HashSet<INode>> CornerNodesOfSubdomains { get; private set; }
        public IReadOnlyDictionary<int, ILinearSystem> LinearSystems { get; }
        public SolverLogger Logger { get; } = new SolverLogger(name);
        public string Name => name;

        public Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider elementMatrixProvider)
        {
            HandleMatrixWillBeSet(); //TODO: temporary solution to avoid this getting called once for each linear system/observable

            var watch = new Stopwatch();
            watch.Start();
            var matrices = new Dictionary<int, IMatrix>();
            foreach (ISubdomain subdomain in model.Subdomains) //TODO: this must be done in parallel
            {
                int s = subdomain.ID;
                IMatrix Kff;
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine($"{this.GetType().Name}: Assembling the free-free stiffness matrix of subdomain {s}");
                    Kff = matrixManagers[s].BuildGlobalMatrix(subdomain.FreeDofOrdering,
                        subdomain.Elements, elementMatrixProvider);
                    linearSystems[s].Matrix = Kff; //TODO: This should be done by the solver not the analyzer. This method should return void.
                }
                else
                {
                    Kff = (IMatrix)(linearSystems[s].Matrix); //TODO: remove the cast
                }
                matrices[s] = Kff;
            }
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);

            this.Initialize(); //TODO: Should this be called by the analyzer? Probably not, since it must be called before DistributeBoundaryLoads().
            return matrices;
        }

        //TODO: There is a lot of repetition between this method and BuildGlobalMatrices
        public Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr)> BuildGlobalSubmatrices(IElementMatrixProvider elementMatrixProvider)
        {
            HandleMatrixWillBeSet(); //TODO: temporary solution to avoid this getting called once for each linear system/observable

            var watch = new Stopwatch();
            watch.Start();
            var matrices = new Dictionary<int, (IMatrix Aff, IMatrixView Afc, IMatrixView Acf, IMatrixView Acc)>();
            foreach (ISubdomain subdomain in model.Subdomains) //TODO: this must be done in parallel
            {
                int s = subdomain.ID;
                if (!subdomain.StiffnessModified)
                {
                    throw new NotImplementedException("This optimization is not implemented");
                }
                if (subdomain.ConstrainedDofOrdering == null)
                {
                    throw new InvalidOperationException("In order to build the matrices corresponding to constrained dofs of,"
                        + $" subdomain {s}, they must have been ordered first.");
                }
                (IMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) =
                    matrixManagers[s].BuildGlobalSubmatrices(subdomain.FreeDofOrdering, subdomain.ConstrainedDofOrdering, 
                    subdomain.Elements, elementMatrixProvider);
                matrices[s] = (Kff, Kfc, Kcf, Kcc);
                linearSystems[s].Matrix = Kff; //TODO: This should be done by the solver not the analyzer. This method should return void.
            }
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);

            this.Initialize(); //TODO: Should this be called by the analyzer? Probably not, since it must be called before DistributeBoundaryLoads().
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
            foreach (ISubdomain subdomain in subdomains.Values)
            {
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine($"{this.GetType().Name}: Clearing saved matrices of subdomain {subdomain.ID}.");
                    matrixManagers[subdomain.ID].Clear();
                }
            }
            flexibility = null;
            preconditioner = null;
            coarseProblemSolver.ClearCoarseProblemMatrix();

            //stiffnessDistribution = null; //WARNING: do not dispose of this. It is updated when BuildGlobalMatrix() is called.
        }

        public void Initialize()
        {
            var watch = new Stopwatch();
            watch.Start();

            // Identify corner nodes
            CornerNodesOfSubdomains = cornerNodeSelection.SelectCornerNodesOfSubdomains(); //TODO: Could this cause change in connectivity?

            // Define boundary / internal dofs
            dofSeparator.DefineGlobalBoundaryDofs(model, CornerNodesOfSubdomains);
            dofSeparator.DefineGlobalCornerDofs(model, CornerNodesOfSubdomains);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int s = subdomain.ID;
                HashSet<INode> cornerNodes = CornerNodesOfSubdomains[s];
                if (subdomain.ConnectivityModified)
                {
                    //TODO: should I cache this somewhere?
                    //TODO: should I use subdomain.Nodes.Except(cornerNodes) instead?
                    IEnumerable<INode> remainderAndConstrainedNodes = subdomain.Nodes.Where(node => !cornerNodes.Contains(node));

                    Debug.WriteLine($"{this.GetType().Name}: Separating and ordering corner-remainder dofs of subdomain {s}");
                    dofSeparator.SeparateCornerRemainderDofs(subdomain, cornerNodes, remainderAndConstrainedNodes);

                    Debug.WriteLine($"{this.GetType().Name}: Reordering internal dofs of subdomain {s}.");
                    matrixManagers[s].ReorderRemainderDofs(dofSeparator, subdomain);

                    Debug.WriteLine($"{this.GetType().Name}: Separating and ordering boundary-internal dofs of subdomain {s}");
                    dofSeparator.SeparateBoundaryInternalDofs(subdomain, remainderAndConstrainedNodes);
                }
            }

            // This must be called after determining the corner dofs of each subdomain
            //TODO: Perhaps the corner dofs of each subdomain should be handled in dofSeparator.DefineGlobalCornerDofs();
            coarseProblemSolver.ReorderCornerDofs(dofSeparator);
            dofSeparator.CalcCornerMappingMatrices(model, CornerNodesOfSubdomains);

            //TODO: B matrices could also be reused in some cases
            // Define lagrange multipliers and boolean matrices. 
            this.lagrangeEnumerator = new FetiDPLagrangeMultipliersEnumerator(crosspointStrategy, dofSeparator);
            if (problemIsHomogeneous) lagrangeEnumerator.DefineBooleanMatrices(model); // optimization in this case
            else lagrangeEnumerator.DefineLagrangesAndBooleanMatrices(model);

            // Log dof statistics
            watch.Stop();
            Logger.LogTaskDuration("Dof ordering", watch.ElapsedMilliseconds);
            int numExpandedDomainFreeDofs = 0;
            foreach (var subdomain in model.Subdomains)
            {
                numExpandedDomainFreeDofs += subdomain.FreeDofOrdering.NumFreeDofs;
            }
            Logger.LogNumDofs("Expanded domain dofs", numExpandedDomainFreeDofs);
            Logger.LogNumDofs("Lagrange multipliers", lagrangeEnumerator.NumLagrangeMultipliers);
            Logger.LogNumDofs("Corner dofs", dofSeparator.NumGlobalCornerDofs);

            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            var Kff = new Dictionary<int, IMatrixView>();
            foreach (int s in linearSystems.Keys) Kff[s] = linearSystems[s].Matrix;
            stiffnessDistribution.Update(Kff);
            subdomainGlobalMapping = new FetiDPSubdomainGlobalMapping(model, dofSeparator, stiffnessDistribution);
        }

        public Dictionary<int, Matrix> InverseSystemMatrixTimesOtherMatrix(Dictionary<int, IMatrixView> otherMatrix)
        {
            throw new NotImplementedException();
        }

        public void OrderDofs(bool alsoOrderConstrainedDofs) //TODO: Only order dofs of subdomains that are modified
        {
            var watch = new Stopwatch();
            watch.Start();

            // Order dofs
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderFreeDofs(model);
            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                if (!subdomain.ConnectivityModified) continue; //TODO: Not sure about this

                matrixManagers[subdomain.ID].HandleDofOrderingWillBeModified();
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
                if (alsoOrderConstrainedDofs) subdomain.ConstrainedDofOrdering = dofOrderer.OrderConstrainedDofs(subdomain);

                // The next must done by the analyzer, so that subdomain.Forces is retained when doing back to back analyses.
                //subdomain.Forces = linearSystem.CreateZeroVector();
            }

            // Log dof statistics
            watch.Stop();
            Logger.LogTaskDuration("Dof ordering", watch.ElapsedMilliseconds);
            Logger.LogNumDofs("Global dofs", globalOrdering.NumGlobalFreeDofs);
        }

        public void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        public void Solve()
        {
            var watch = new Stopwatch();
            foreach (var linearSystem in linearSystems.Values)
            {
                if (linearSystem.SolutionConcrete == null) linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
            }

            // Separate the force vector
            watch.Start();
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
            watch.Stop();
            Logger.LogTaskDuration("Separating vectors & matrices", watch.ElapsedMilliseconds);
            watch.Reset();

            if (isStiffnessModified)
            {
                // Separate the stiffness matrix
                watch.Start();
                foreach (int s in subdomains.Keys)
                {
                    if (!subdomains[s].StiffnessModified) continue;
                    IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];
                    int[] remainderDofs = dofSeparator.RemainderDofIndices[s];
                    int[] cornerDofs = dofSeparator.CornerDofIndices[s];
                    matrices.ExtractKrr(remainderDofs);
                    matrices.ExtractKcrKrc(cornerDofs, remainderDofs);
                    matrices.ExtractKcc(cornerDofs);
                }
                watch.Stop();
                Logger.LogTaskDuration("Separating vectors & matrices", watch.ElapsedMilliseconds);


                // Reorder internal dofs if needed by the preconditioner. TODO: Should I have done this previously in Initialize()?
                watch.Start();
                if (preconditionerFactory.ReorderInternalDofsForFactorization)
                {
                    foreach (ISubdomain subdomain in model.Subdomains)
                    {
                        if (subdomain.ConnectivityModified)
                        {
                            Debug.WriteLine($"{this.GetType().Name}: Reordering internal dofs of subdomain {subdomain.ID}.");
                            matrixManagers[subdomain.ID].ReorderInternalDofs(dofSeparator, subdomain);
                        }
                    }
                }

                // Calculate the preconditioner before factorizing each subdomain's Kff 
                preconditioner = preconditionerFactory.CreatePreconditioner(model, stiffnessDistribution, dofSeparator,
                    lagrangeEnumerator, matrixManagersGeneral);
                watch.Stop();
                Logger.LogTaskDuration("Calculating preconditioner", watch.ElapsedMilliseconds);

                // Factorize each subdomain's Krr
                watch.Restart();
                foreach (int s in subdomains.Keys)
                {
                    if (!subdomains[s].StiffnessModified) continue;
                    //TODO: If I can reuse Krr, I can also reuse its factorization. Therefore this must be inPlace. In contrast, FETI-1 needs Kff intact for Stiffness distribution, in the current design).
                    Debug.WriteLine(
                        $"{this.GetType().Name}: Inverting the remainder-remainder stiffness matrix of subdomain {s} in place.");
                    matrixManagers[s].InvertKrr(true);
                }
                watch.Stop();
                Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);

                // Define FETI-DP flexibility matrices
                watch.Restart();
                flexibility = new FetiDPFlexibilityMatrix(dofSeparator, lagrangeEnumerator, matrixManagers);

                // Static condensation of remainder dofs (Schur complement).
                coarseProblemSolver.CreateAndInvertCoarseProblemMatrix(CornerNodesOfSubdomains, dofSeparator, matrixManagers);
                watch.Stop();
                Logger.LogTaskDuration("Setting up interface problem", watch.ElapsedMilliseconds);
                watch.Reset();

                isStiffnessModified = false;
            }

            // Static condensation for the force vectors
            watch.Start();
            Vector globalFcStar = coarseProblemSolver.CreateCoarseProblemRhs(dofSeparator, matrixManagers, fr, fbc);

            // Calculate the rhs vectors of the interface system
            Vector dr = CalcDisconnectedDisplacements(fr);
            double globalForcesNorm = CalcGlobalForcesNorm();
            watch.Stop();
            Logger.LogTaskDuration("Setting up interface problem", watch.ElapsedMilliseconds);

            // Solve the interface problem
            watch.Restart();
            (Vector lagranges, Vector uc) = interfaceProblemSolver.SolveInterfaceProblem(flexibility, preconditioner, 
                coarseProblemSolver, globalFcStar, dr, globalForcesNorm, Logger);
            watch.Stop();
            Logger.LogTaskDuration("Solving interface problem", watch.ElapsedMilliseconds);

            // Calculate the displacements of each subdomain
            watch.Restart();
            Dictionary<int, Vector> actualDisplacements = CalcActualDisplacements(lagranges, uc,  fr);
            foreach (var idSystem in linearSystems) idSystem.Value.SolutionConcrete = actualDisplacements[idSystem.Key];
            watch.Stop();
            Logger.LogTaskDuration("Calculate displacements from lagrange multipliers", watch.ElapsedMilliseconds);

            Logger.IncrementAnalysisStep();
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

        public class Builder
        {
            private readonly ICornerNodeSelection cornerNodeSelection;
            private readonly IFetiDPSubdomainMatrixManagerFactory matrixManagerFactory; 

            public Builder(ICornerNodeSelection cornerNodeSelection, 
                IFetiDPSubdomainMatrixManagerFactory matrixManagerFactory)
            {
                this.cornerNodeSelection = cornerNodeSelection;
                this.matrixManagerFactory = matrixManagerFactory;
            }

            //TODO: We need to specify the ordering for remainder and possibly internal dofs, while IDofOrderer only works for free dofs.
            public IDofOrderer DofOrderer { get; set; } =
                new ReusingDofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public IFetiDPInterfaceProblemSolver InterfaceProblemSolver { get; set; } 
                = new FetiDPInterfaceProblemSolver.Builder().Build();
            public IFetiPreconditionerFactory PreconditionerFactory { get; set; } = new LumpedPreconditioner.Factory();
            public bool ProblemIsHomogeneous { get; set; } = true;

            public FetiDPSolver BuildSolver(IStructuralModel model)
                => new FetiDPSolver(model, cornerNodeSelection, matrixManagerFactory, DofOrderer, PreconditionerFactory,
                        ProblemIsHomogeneous, InterfaceProblemSolver);            
        }
    }
}
