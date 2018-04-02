using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Logging
{
    //TODO: add MatrixMarket, etc classes
    public class SparseMatrixFormatting
    {
        public SparseMatrixFormatting()
        {
        }

        public string FormatNonZeroEntry(int row, int col, double value)
        {
            return $"A[{row}, {col}] = {value}\n";
        }
    }
}
