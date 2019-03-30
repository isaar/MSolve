using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition
{
    public interface IDofSeparator
    {
        Dictionary<INode, DOFType[]> GlobalBoundaryDofs { get; }
    }
}
