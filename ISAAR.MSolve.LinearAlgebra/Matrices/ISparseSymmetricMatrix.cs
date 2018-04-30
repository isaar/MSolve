using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public interface ISparseSymmetricMatrix: ISymmetricMatrix, ISparseMatrix
    {
        int CountNonZerosSuperDiagonal();
        IEnumerable<(int row, int col, double value)> EnumerateNonZerosSuperDiagonal(); //Perhaps the matrix itself should be the target of foreach
    }
}
