using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.XFEM.Utilities
{
    static class MatrixUtilities
    {
        public static void AddPartialToTotalMatrix(Matrix2D<double> partialMatrix, Matrix2D<double> totalMatrix)
        {
            Debug.Assert(partialMatrix.Rows == totalMatrix.Rows, "Non matching rows");
            Debug.Assert(partialMatrix.Columns == totalMatrix.Columns, "Non matching columns");
            for (int row = 0; row < totalMatrix.Rows; ++row)
            {
                for (int col = 0; col < totalMatrix.Columns; ++col)
                {
                    totalMatrix[row, col] += partialMatrix[row, col];
                }
            }
        }

        public static void AddPartialToSymmetricTotalMatrix(Matrix2D<double> partialMatrix,
            SymmetricMatrix2D<double> totalMatrix)
        {
            Debug.Assert(partialMatrix.Rows == totalMatrix.Rows, "Non matching rows");
            Debug.Assert(partialMatrix.Columns == totalMatrix.Columns, "Non matching columns");
            for (int row = 0; row < totalMatrix.Rows; ++row)
            {
                for (int col = row; col < totalMatrix.Columns; ++col)
                {
                    totalMatrix[row, col] += partialMatrix[row, col];
                }
            }
        }
    }
}
