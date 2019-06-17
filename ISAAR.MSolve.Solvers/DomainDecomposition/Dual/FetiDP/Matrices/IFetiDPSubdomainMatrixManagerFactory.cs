using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices
{
    public interface IFetiDPSubdomainMatrixManagerFactory
    {
        //TODO: This is in a different namespace
        IFetiDPCoarseProblemSolver CreateCoarseProblemSolver(IReadOnlyList<ISubdomain> subdomains); 

        IFetiDPSubdomainMatrixManager CreateMatricesManager(ISubdomain subdomain);
    }
}
