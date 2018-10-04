using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Commons
{
    public interface ILinearSystem_v2
    {
        //TODO: this is error prone. The implementation should manage the state, by funneling access to the matrix.
        bool IsMatrixModified { get; set; }

        IMatrix Matrix { get; }
        IVector RhsVector { get; }
        IVector Solution { get; }
    }
}
