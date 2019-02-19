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
        private Dictionary<int, IMatrix> capacityMatrices, conductivityMatrices;
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

        private IDictionary<int, IMatrix> CapacityMatrices
        {
            get
            {
                if (capacityMatrices == null) BuildCapacityMatrices();
                return capacityMatrices;
            }
        }

        private IDictionary<int, IMatrix> ConductivityMatrices
        {
            get
            {
                if (conductivityMatrices == null) BuildConductivityMatrices();
                //else RebuildConductivityMatrices();
                return conductivityMatrices;
            }
        }

        private void BuildConductivityMatrices()
        {
            conductivityMatrices = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                conductivityMatrices.Add(subdomain.ID, solver.BuildGlobalMatrix(subdomain, conductivityProvider));
            }
        }

        private void RebuildConductivityMatrices()
        {
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                if (subdomain.MaterialsModified)
                {
                    conductivityMatrices[subdomain.ID] = solver.BuildGlobalMatrix(subdomain, conductivityProvider);
                    subdomain.ResetMaterialsModifiedProperty();
                }
            }
        }

        private void BuildCapacityMatrices()
        {
            capacityMatrices = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                capacityMatrices.Add(subdomain.ID, solver.BuildGlobalMatrix(subdomain, capacityProvider));
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

            conductivityMatrices = null;
            capacityMatrices = null;
        }
        #endregion 

        #region IImplicitIntegrationProvider Members

        public IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients,
            ISubdomain_v2 subdomain)
        {
            // The effective matrix should not overwrite the conductivity matrix. 
            // In a dynamic analysis that is not purely implicit we need the conductivity matrix.
            int id = subdomain.ID;
            return ConductivityMatrices[id].LinearCombination(coefficients.Stiffness, CapacityMatrices[id], coefficients.Mass);
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

            model.AssignNodalLoads(); // Time-independent nodal loads
            model.AssignTimeDependentNodalLoads(timeStep); // Time-dependent nodal loads

            var rhsVectors = new Dictionary<int, IVector>();
            foreach (ISubdomain_v2 subdomain in model.Subdomains) rhsVectors.Add(subdomain.ID, subdomain.Forces.Copy());
            return rhsVectors;
        }

        public IVector MassMatrixVectorProduct(ISubdomain_v2 subdomain, IVectorView vector)
            => this.CapacityMatrices[subdomain.ID].Multiply(vector);

        //TODO: Ok this is weird. These methods should be named Second/First/ZeroOrderCoefficientTimesVector()
        public IVector DampingMatrixVectorProduct(ISubdomain_v2 subdomain, IVectorView vector)
            => this.ConductivityMatrices[subdomain.ID].Multiply(vector);

        #endregion

        #region IStaticProvider Members

        public IMatrixView CalculateMatrix(ISubdomain_v2 subdomain)
        {
            if (conductivityMatrices == null) BuildConductivityMatrices();
            return conductivityMatrices[subdomain.ID];
        }

        #endregion

        #region INonLinearProvider Members

        public double CalculateRhsNorm(IVectorView rhs) => rhs.Norm2();

        public void ProcessInternalRhs(ISubdomain_v2 subdomain, IVectorView solution, IVector rhs) { }

        #endregion
    }
}
