using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition
{
    public interface IDofSeparator
    {
        Dictionary<INode, IDofType[]> GlobalBoundaryDofs { get; }
    }
}
