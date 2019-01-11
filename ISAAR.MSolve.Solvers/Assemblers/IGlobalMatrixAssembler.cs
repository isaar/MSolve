using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: not sure this interface is required
namespace ISAAR.MSolve.Solvers.Assemblers
{
    public interface IGlobalMatrixAssembler<TMatrix>
        where TMatrix : IMatrix
    {
        TMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement_v2> elements,
            IElementMatrixProvider_v2 matrixProvider);

        /// <summary>
        /// It will be called before modifying the current freedom degree ordering.
        /// </summary>
        void HandleDofOrderingWillBeModified();
    }
}
