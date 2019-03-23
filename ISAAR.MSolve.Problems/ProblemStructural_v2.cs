using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: Usually the LinearSystem is passed in, but for GetRHSFromHistoryLoad() it is stored as a field. Decide on one method.
//TODO: I am not too fond of the provider storing global sized matrices.
namespace ISAAR.MSolve.Problems
{
    public class ProblemStructural_v2 : IImplicitIntegrationProvider_v2, IStaticProvider_v2, INonLinearProvider_v2
    {
        private Dictionary<int, IMatrix> ms, cs, ks;
        private readonly IStructuralModel_v2 model;
        private readonly ISolver_v2 solver;
        private IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        private ElementStructuralStiffnessProvider_v2 stiffnessProvider = new ElementStructuralStiffnessProvider_v2();
        private ElementStructuralMassProvider_v2 massProvider = new ElementStructuralMassProvider_v2();
        private ElementStructuralDampingProvider_v2 dampingProvider = new ElementStructuralDampingProvider_v2();

        public ProblemStructural_v2(IStructuralModel_v2 model, ISolver_v2 solver)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.DirichletLoadsAssembler = new DirichletEquivalentLoadsStructural(stiffnessProvider);
        }

        //public double AboserberE { get; set; }
        //public double Aboseberv { get; set; }

        public IDirichletEquivalentLoadsAssembler DirichletLoadsAssembler { get; } 

        private IDictionary<int, IMatrix> Ms
        {
            get
            {
                if (ms == null) BuildMs();
                return ms;
            }
        }

        private IDictionary<int, IMatrix> Cs
        {
            get
            {
                if (cs == null) BuildCs();
                return cs;
            }
        }

        private IDictionary<int, IMatrix> Ks
        {
            get
            {
                if (ks == null)
                    BuildKs();
                else
                {
                    //TODO I am not too fond of side effects, especially in getters
                    RebuildKs(); // This is the same but also resets the material modified properties. 
                }
                return ks;
            }
        }

        private void BuildKs()
        {
            ks = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                ks.Add(subdomain.ID, solver.BuildGlobalMatrix(subdomain, stiffnessProvider));
            }
        }

        private void RebuildKs()
        {
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                if (subdomain.MaterialsModified)
                {
                    ks[subdomain.ID] = solver.BuildGlobalMatrix(subdomain, stiffnessProvider);
                    subdomain.ResetMaterialsModifiedProperty();
                }
            }
        }

        private void BuildMs()
        {
            ms = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                ms.Add(subdomain.ID, solver.BuildGlobalMatrix(subdomain, massProvider));
            }
        }

        //TODO: With Rayleigh damping, C is more efficiently built using linear combinations of global K, M, 
        //      instead of building and assembling element k, m matrices.
        private void BuildCs()
        {
            cs = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                cs.Add(subdomain.ID, solver.BuildGlobalMatrix(subdomain, dampingProvider));
            }
        }

        #region IAnalyzerProvider Members
        public void Reset()
        {
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                foreach (IElement_v2 element in subdomain.Elements)
                {
                    ((IFiniteElement_v2)element.ElementType).ClearMaterialState();
                }
            }

            cs = null;
            ks = null;
            ms = null;
        }
        #endregion 

        #region IImplicitIntegrationProvider Members

        public IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients, 
            ISubdomain_v2 subdomain)
        {
            //TODO: 1) Why do we want Ks to be built only if it has not been factorized? 
            //      2) When calling Ks[id], the matrix will be built anyway, due to the annoying side effects of the property.
            //         Therefore, if the matrix was indeed factorized it would be built twice!
            //      3) The provider should be decoupled from solver logic, such as knowing if the matrix is factorized. Knowledge
            //         that the matrix has been altered by the solver could be implemented by observers, if necessary.
            //      4) The analyzer should decide when global matrices need to be rebuilt, not the provider.
            //      5) The need to rebuild the system matrix if the solver has modified it might be avoidable if the analyzer 
            //         uses the appropriate order of operations. However, that may not always be possible. Such a feature 
            //         (rebuild or store) is nice to have. Whow would be responsible, the solver, provider or assembler?
            //      6) If the analyzer needs the system matrix, then it can call solver.PreventFromOverwritingMatrix(). E.g.
            //          explicit dynamic analyzers would need to do that.
            //if (linearSystem.IsMatrixOverwrittenBySolver) BuildKs();

            int id = subdomain.ID;
            IMatrix matrix = this.Ks[id];
            matrix.LinearCombinationIntoThis(coefficients.Stiffness, Ms[id], coefficients.Mass);
            matrix.AxpyIntoThis(Cs[id], coefficients.Damping);
            return this.Ks[id];
        }

        public void ProcessRhs(ImplicitIntegrationCoefficients coefficients, ISubdomain_v2 subdomain, IVector rhs)
        {
            // Method intentionally left empty.
        }

        public IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep)
        {
            var d = new Dictionary<int, IVector>();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                d.Add(linearSystem.Subdomain.ID, linearSystem.CreateZeroVector());
            }

            if (model.MassAccelerationHistoryLoads.Count > 0)
            {
                List<MassAccelerationLoad> m = new List<MassAccelerationLoad>(model.MassAccelerationHistoryLoads.Count);
                foreach (IMassAccelerationHistoryLoad l in model.MassAccelerationHistoryLoads)
                    m.Add(new MassAccelerationLoad() { Amount = l[timeStep], DOF = l.DOF });

                foreach (ISubdomain_v2 subdomain in model.Subdomains)
                {
                    int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                    foreach ((INode node, DOFType dofType, int subdomainDofIdx) in subdomain.DofOrdering.FreeDofs)
                    {
                        int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                        foreach (var l in m)
                        {
                            if (dofType == l.DOF) d[subdomain.ID].Set(globalDofIdx, l.Amount);
                        }
                    }

                    //foreach (var nodeInfo in subdomain.GlobalNodalDOFsDictionary)
                    //{
                    //    foreach (var dofPair in nodeInfo.Value)
                    //    {
                    //        foreach (var l in m)
                    //        {
                    //            if (dofPair.Key == l.DOF && dofPair.Value != -1)
                    //            {
                    //                d[subdomain.ID].Set(dofPair.Value, l.Amount);
                    //            }
                    //        }
                    //    }
                    //}
                }
            }

            //foreach (ElementMassAccelerationHistoryLoad load in model.ElementMassAccelerationHistoryLoads)
            //{
            //    MassAccelerationLoad hl = new MassAccelerationLoad() { Amount = load.HistoryLoad[timeStep] * 564000000, DOF = load.HistoryLoad.DOF };
            //    load.Element.Subdomain.AddLocalVectorToGlobal(load.Element,
            //        load.Element.ElementType.CalculateAccelerationForces(load.Element, (new MassAccelerationLoad[] { hl }).ToList()),
            //        load.Element.Subdomain.Forces);
            //}

            return d;
        }

        public IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep)
        {
            var d = new Dictionary<int, IVector>();

            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                d.Add(linearSystem.Subdomain.ID, linearSystem.CreateZeroVector());
            }

            return d;
        }

        public IDictionary<int, IVector> GetRhsFromHistoryLoad(int timeStep)
        {
            foreach (ISubdomain_v2 subdomain in model.Subdomains) subdomain.Forces.Clear(); //TODO: this is also done by model.AssignLoads()

            model.AssignLoads();
            model.AssignMassAccelerationHistoryLoads(timeStep);

            var rhsVectors = new Dictionary<int, IVector>();
            foreach (Subdomain_v2 subdomain in model.Subdomains) rhsVectors.Add(subdomain.ID, subdomain.Forces);
            return rhsVectors;
        }

        public IVector MassMatrixVectorProduct(ISubdomain_v2 subdomain, IVectorView vector)
            => this.Ms[subdomain.ID].Multiply(vector);

        public IVector DampingMatrixVectorProduct(ISubdomain_v2 subdomain, IVectorView vector)
            => this.Cs[subdomain.ID].Multiply(vector);

        #endregion

        #region IStaticProvider Members

        public IMatrixView CalculateMatrix(ISubdomain_v2 subdomain)
        {
            if (ks == null) BuildKs();
            return ks[subdomain.ID];
        }

        #endregion

        #region INonLinearProvider Members

        public double CalculateRhsNorm(IVectorView rhs) => rhs.Norm2();

        public void ProcessInternalRhs(ISubdomain_v2 subdomain, IVectorView solution, IVector rhs) {}

        #endregion
    }
}
