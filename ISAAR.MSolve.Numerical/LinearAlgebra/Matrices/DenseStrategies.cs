using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    internal static class DenseStrategies
    {
        internal static bool AreEqual(IIndexable2D matrix1, IIndexable2D matrix2, double tolerance = 1e-13)
        {
            if ((matrix1.NumRows != matrix2.NumRows) || (matrix1.NumColumns != matrix2.NumColumns)) return false;
            var comparer = new ValueComparer(1e-13);
            for (int j = 0; j < matrix1.NumColumns; ++j)
            {
                for (int i = 0; i < matrix1.NumRows; ++i)
                {
                    if (!comparer.AreEqual(matrix1[i, j], matrix2[i, j])) return false;
                }
            }
            return true;
        }

        internal static Matrix DoEntrywise(IIndexable2D matrix1, IIndexable2D matrix2, 
            Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(matrix1, matrix2);
            var result = Matrix.CreateZero(matrix1.NumRows, matrix1.NumColumns);
            for (int j = 0; j < matrix1.NumColumns; ++j)
            {
                for (int i = 0; i < matrix1.NumRows; ++i)
                {
                    result[i, j] = binaryOperation(matrix1[i, j], matrix2[i, j]);
                }
            }
            return result;
        }

        internal static Matrix Transpose(ITransposable matrix)
        {
            var result = Matrix.CreateZero(matrix.NumColumns, matrix.NumRows);
            for (int j = 0; j < matrix.NumColumns; ++j)
            {
                for (int i = 0; i < matrix.NumRows; ++i)
                {
                    result[j, i] = matrix[i, j];
                }
            }
            return result;
        }
    }
}
