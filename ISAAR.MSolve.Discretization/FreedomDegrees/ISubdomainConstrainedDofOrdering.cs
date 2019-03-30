using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: many of these operations would be more efficient if done simultaneously for free and constrained dofs.
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public interface ISubdomainConstrainedDofOrdering
    {
        DofTable ConstrainedDofs { get; }
        int NumConstrainedDofs { get; }
        (int[] elementDofIndices, int[] subdomainDofIndices) MapConstrainedDofsElementToSubdomain(IElement_v2 element);
    }
}
