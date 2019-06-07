using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices
{
    public interface IFetiDPSubdomainMatrixManagerFactory
    {
        IFetiDPSubdomainMatrixManager CreateMatricesManager(ISubdomain subdomain);
    }
}
