using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    public interface ISparsityPattern
    {
        int NumColumns { get; }
        int NumRows { get; }
        int CountNonZeros();
        IEnumerable<(int row, int col)> EnumerateNonZeros();
        bool IsNonZero(int rowIdx, int colIdx);
    }
}
