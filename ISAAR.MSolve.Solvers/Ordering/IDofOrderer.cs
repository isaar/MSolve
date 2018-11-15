using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.Ordering
{
    public interface IDofOrderer
    {
        bool AreDofsOrdered { get; } //TODO: I do not like this. It would be better to use a factory instead.
        DofTable FreeDofs { get; }
        int NumFreeDofs { get; }

        IVector ExtractVectorElementFromGlobal(IElement element, IVectorView globalVector);
        IReadOnlyDictionary<int, int> MapFreeDofsElementToGlobal(IElement element);
        void OrderDofs(ISubdomain subdomain);
    }
}
