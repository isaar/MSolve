using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Assemblers
{
    public interface IGlobalMatrixAssembler
    {
        void BuildGlobalMatrix(IEnumerable<IElement> elements, FreeDofOrderer dofOrderer, IElementMatrixProvider matrixProvider);
    }
}
