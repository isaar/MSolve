using System;
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
    public class ProblemThermal_v2 : IImplicitIntegrationProvider_v2, IStaticProvider_v2, INonLinearProvider_v2
    {
        private Dictionary<int, IMatrix> capacity, conductivityFreeFree;
        private Dictionary<int, IMatrixView> conductivityFreeConstr, conductivityConstrFree, conductivityConstrConstr;
        private readonly Model_v2 model;
        private readonly ISolver_v2 solver;
        private IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        private ElementStructuralStiffnessProvider_v2 conductivityProvider = new ElementStructuralStiffnessProvider_v2();
        private ElementStructuralMassProvider_v2 capacityProvider = new ElementStructuralMassProvider_v2();

        public ProblemThermal_v2(Model_v2 model, ISolver_v2 solver)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.DirichletLoadsAssembler = new DirichletEquivalentLoadsStructural(conductivityProvider);
        }

        public IDirichletEquivalentLoadsAssembler DirichletLoadsAssembler { get; }

        private IDictionary<int, IMatrix> Capacity
        {
            get
            {
                if (capacity == null) BuildCapacity();
                return capacity;
            }
        }

        private IDictionary<int, IMatrix> Conductivity
        {
            get
            {
                if (conductivityFreeFree == null) BuildConductivityFreeFree();
                //else RebuildConductivityMatrices();
                return conductivityFreeFree;
            }
        }

        private void BuildConductivityFreeFree() => conductivityFreeFree = solver.BuildGlobalMatrices(conductivityProvider);

        private void BuildConductivitySubmatrices()
        {
            Dictionary<int, (IMatrix Cff, IMatrixView Cfc, IMatrixView Ccf, IMatrixView Ccc)> matrices =
                solver.BuildGlobalSubmatrices(conductivityProvider);

            conductivityFreeFree = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            conductivityFreeConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            conductivityConstrFree = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            conductivityConstrConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                conductivityFreeFree.Add(id, matrices[id].Cff);
                conductivityFreeConstr.Add(id, matrices[id].Cfc);
                conductivityConstrFree.Add(id, matrices[id].Ccf);
                conductivityConstrConstr.Add(id, matrices[id].Ccc);
            }
        }

        private void RebuildConductivityFreeFree()
        {
            //TODO: This will rebuild all the stiffnesses of all subdomains, if even one subdomain has MaterialsModified = true.
            //      Optimize this, by passing a flag foreach subdomain to solver.BuildGlobalSubmatrices().

            bool mustRebuild = false;
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                if (subdomain.MaterialsModified)
                {
                    mustRebuild = true;
                    break;
                }
            }
            if (mustRebuild) conductivityFreeFree = solver.BuildGlobalMatrices(conductivityProvider);
            foreach (ISubdomain_v2 subdomain in model.Subdomains) subdomain.ResetMaterialsModifiedProperty();
        }

        private void BuildCapacity() => capacity = solver.BuildGlobalMatrices(capacityProvider);

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

            conductivityFreeFree = null;
            conductivityConstrFree = null;
            conductivityConstrConstr = null;
            capacity = null;
        }
        #endregion 

        #region IImplicitIntegrationProvider Members

        public IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients,
            ISubdomain_v2 subdomain)
        {
            // The effective matrix should not overwrite the conductivity matrix. 
            // In a dynamic analysis that is not purely implicit we need the conductivity matrix.
            int id = subdomain.ID;
            return Conductivity[id].LinearCombination(coefficients.Stiffness, Capacity[id], coefficients.Mass);
        }

        public void ProcessRhs(ImplicitIntegrationCoefficients coefficients, ISubdomain_v2 subdomain, IVector rhs)
        {
            // Method intentionally left empty.
        }

        public IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep)
        {
            throw new InvalidOperationException("This is does not make sense in explicit methods for first order equations");
        }

        public IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep)
        {
            throw new InvalidOperationException("This is not needed in explicit methods for first order equations");
        }

        public IDictionary<int, IVector> GetRhsFromHistoryLoad(int timeStep)
        {
            foreach (ISubdomain_v2 subdomain in model.Subdomains) subdomain.Forces.Clear(); //TODO: this is also done by model.AssignLoads()

            model.AssignNodalLoads(solver.DistributeNodalLoads); // Time-independent nodal loads
            model.AssignTimeDependentNodalLoads(timeStep, solver.DistributeNodalLoads); // Time-dependent nodal loads

            var rhsVectors = new Dictionary<int, IVector>();
            foreach (ISubdomain_v2 subdomain in model.Subdomains) rhsVectors.Add(subdomain.ID, subdomain.Forces.Copy());
            return rhsVectors;
        }

        public IVector MassMatrixVectorProduct(ISubdomain_v2 subdomain, IVectorView vector)
            => this.Capacity[subdomain.ID].Multiply(vector);

        //TODO: Ok this is weird. These methods should be named Second/First/ZeroOrderCoefficientTimesVector()
        public IVector DampingMatrixVectorProduct(ISubdomain_v2 subdomain, IVectorView vector)
            => this.Conductivity[subdomain.ID].Multiply(vector);

        #endregion

        #region IStaticProvider Members

        public IMatrixView CalculateMatrix(ISubdomain_v2 subdomain)
        {
            if (conductivityFreeFree == null) BuildConductivityFreeFree();
            return conductivityFreeFree[subdomain.ID];
        }

        public (IMatrixView matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr) CalculateSubMatrices(ISubdomain_v2 subdomain)
        {
            int id = subdomain.ID;
            if ((conductivityFreeFree == null) || (conductivityFreeConstr == null) 
                || (conductivityConstrFree == null) || (conductivityConstrConstr == null))
            {
                BuildConductivitySubmatrices();
            }
            return (conductivityFreeFree[id], conductivityFreeConstr[id], 
                conductivityConstrFree[id], conductivityConstrConstr[id]);
        }
        #endregion

        #region INonLinearProvider Members

        public double CalculateRhsNorm(IVectorView rhs) => rhs.Norm2();

        public void ProcessInternalRhs(ISubdomain_v2 subdomain, IVectorView solution, IVector rhs) { }

        #endregion
    }
}
