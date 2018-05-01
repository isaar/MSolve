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

        //Perhaps the matrix itself should be the target of foreach. Is this lazy evaluation? Do I need it?
        IEnumerable<(int row, int col, double value)> EnumerateNonZeros(); 

        SparseFormat GetSparseFormat();
    }
}
