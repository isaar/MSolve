using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Output;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public interface ISparseMatrix: IIndexable2D
    {
        int CountNonZeros();
        IEnumerable<(int row, int col, double value)> EnumerateNonZeros();
        SparseFormat GetSparseFormat();
    }
}
