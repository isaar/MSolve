using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public interface ISquareMatrix: IMatrix, ISquareMatrixView
    {
        // This is the common usecase for storing L and U together. ISquareMatrixView should return them separately.
        ILUFactorization FactorLUIntoThis(bool ldu = false); 
        //void FactorLDLIntoThis(); This belongs in SymmetricMatrix
    }
}
