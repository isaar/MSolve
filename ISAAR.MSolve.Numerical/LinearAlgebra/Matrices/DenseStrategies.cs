using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

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

        internal static double[,] CopyToArray2D(IIndexable2D matrix)
        {
            var result = new double[matrix.NumRows, matrix.NumColumns];
            for (int j = 0; j < matrix.NumColumns; ++j)
            {
                for (int i = 0; i < matrix.NumRows; ++i)
                {
                    result[i, j] = matrix[i, j];
                }
            }
            return result;
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

        internal static Matrix Transpose(IMatrixView matrix)
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

        internal static Matrix Multiply(IMatrixView matrix1, IMatrixView matrix2, bool transpose1, bool transpose2)
        {
            if (transpose1)
            {
                if (transpose2)
                {
                    Preconditions.CheckMultiplicationDimensions(matrix1.NumRows, matrix2.NumColumns);
                    var result = Matrix.CreateZero(matrix1.NumColumns, matrix2.NumRows);
                    for (int i = 0; i < result.NumRows; ++i)
                    {
                        for (int j = 0; j < result.NumColumns; ++j)
                        {
                            for (int k = 0; k < matrix1.NumRows; ++k)
                            {
                                result[i, j] = matrix1[k, i] * matrix2[j, k];
                            }
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(matrix1.NumRows, matrix2.NumRows);
                    var result = Matrix.CreateZero(matrix1.NumColumns, matrix2.NumColumns);
                    for (int i = 0; i < result.NumRows; ++i)
                    {
                        for (int j = 0; j < result.NumColumns; ++j)
                        {
                            for (int k = 0; k < matrix1.NumRows; ++k)
                            {
                                result[i, j] = matrix1[k, i] * matrix2[k, j];
                            }
                        }
                    }
                    return result;
                }
            }
            else
            {
                if (transpose2)
                {
                    Preconditions.CheckMultiplicationDimensions(matrix1.NumColumns, matrix2.NumColumns);
                    var result = Matrix.CreateZero(matrix1.NumRows, matrix2.NumRows);
                    for (int i = 0; i < result.NumRows; ++i)
                    {
                        for (int j = 0; j < result.NumColumns; ++j)
                        {
                            for (int k = 0; k < matrix1.NumColumns; ++k)
                            {
                                result[i, j] = matrix1[i, k] * matrix2[j, k];
                            }
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(matrix1.NumColumns, matrix2.NumRows);
                    var result = Matrix.CreateZero(matrix1.NumRows, matrix2.NumColumns);
                    for (int i = 0; i < result.NumRows; ++i)
                    {
                        for (int j = 0; j < result.NumColumns; ++j)
                        {
                            for (int k = 0; k <matrix1.NumColumns; ++k)
                            {
                                result[i, j] = matrix1[i, k] * matrix2[k, j];
                            }
                        }
                    }
                    return result;
                }
            }
        }

        internal static VectorMKL Multiply(IMatrixView matrix, IVectorView vector, bool transposeMatrix)
        {
            if (transposeMatrix)
            {
                Preconditions.CheckMultiplicationDimensions(matrix.NumRows, vector.Length);
                var result = VectorMKL.CreateZero(matrix.NumColumns);
                for (int i = 0; i < result.Length; ++i)
                {
                    for (int j = 0; j < vector.Length; ++j)
                    {
                        result[i] = matrix[j, i] * vector[j];
                    }
                }
                return result;
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, vector.Length);
                var result = VectorMKL.CreateZero(matrix.NumRows);
                for (int i = 0; i < result.Length; ++i)
                {
                    for (int j = 0; j < vector.Length; ++j)
                    {
                        result[i] = matrix[j, i] * vector[j];
                    }
                }
                return result;
            }
        }
    }
}
