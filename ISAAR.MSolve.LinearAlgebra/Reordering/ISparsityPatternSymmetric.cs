using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    public interface ISparsityPatternSymmetric: ISparsityPattern
    {
        int CountNonZerosUpper();
        IEnumerable<(int row, int col)> EnumerateNonZerosUpper();
    }
}
