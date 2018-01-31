using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public interface ISquareMatrixView: IMatrixView
    {
        int Order { get; }

        IVector ExtractDiagonal(); // This is only for square matrices
        double CalcDeterminant();
        ILUFactorization FactorLU(bool ldu = false); // For now L,U are inside the same matrix
        ISquareMatrix Invert();
    }
}
