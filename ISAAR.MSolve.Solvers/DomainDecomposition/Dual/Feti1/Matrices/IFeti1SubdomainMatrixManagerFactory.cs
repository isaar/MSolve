using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices
{
    public interface IFeti1SubdomainMatrixManagerFactory
    {
        IFeti1SubdomainMatrixManager CreateMatricesManager(ISubdomain subdomain);
    }
}
