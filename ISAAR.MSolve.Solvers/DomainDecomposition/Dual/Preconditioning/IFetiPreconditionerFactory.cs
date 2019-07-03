using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Also add ReorderInternalDofsForMultiplication and ReorderBoundaryDofsForMultiplication
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public interface IFetiPreconditionerFactory
    {
        bool ReorderInternalDofsForFactorization { get; }

        IFetiPreconditioner CreatePreconditioner(IStructuralModel model, IStiffnessDistribution stiffnessDistribution,
            IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangeEnumerator, 
            Dictionary<int, IFetiSubdomainMatrixManager> matrixManagers);
    }
}
