using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public interface IDofOrdering
    {
        DofTable FreeDofs { get; }
        int NumFreeDofs { get; }

        IVector ExtractVectorElementFromGlobal(IElement element, IVectorView globalVector);
        IReadOnlyDictionary<int, int> MapFreeDofsElementToGlobal(IElement element);
    }
}
