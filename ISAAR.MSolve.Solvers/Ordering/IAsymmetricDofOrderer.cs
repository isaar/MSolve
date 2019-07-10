using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Solvers.Ordering
{
    public interface IAsymmetricDofOrderer
    {
        ISubdomainConstrainedDofOrdering OrderConstrainedDofs(ISubdomain subdomain);
        IGlobalFreeDofOrdering OrderFreeDofs(IStructuralAsymmetricModel model);
    }
}
