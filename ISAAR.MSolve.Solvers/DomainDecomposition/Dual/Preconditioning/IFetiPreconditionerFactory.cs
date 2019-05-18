using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public interface IFetiPreconditionerFactory
    {
        IFetiPreconditioner CreatePreconditioner(IStiffnessDistribution stiffnessDistribution, IDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangeEnumerator, Dictionary<int, IMatrixView> stiffnessMatrices);
    }
}
