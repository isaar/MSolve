using System;
using ISAAR.MSolve.LinearAlgebra.Exceptions;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Utilities
{
    /// <summary>
    /// Linear algebra operations to be tested against. These are naive implementations for dense matrices and vectors, 
    /// thus they are much more reliable. Still getting the result from another software (e.g. matlab) would be preferable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MatrixOperations
    {
        public static double[] LinearCombination(double scalar1, double[] vector1, double scalar2, double[] vector2)
        {
            int n = vector1.Length;
            if (vector2.Length != n) throw new NonMatchingDimensionsException("Cannot add arrays with different length");
            var c = new double[n];
            for (int i = 0; i < n; ++i)
            {
                c[i] = scalar1 * vector1[i] + scalar2 * vector2[i];
            }
            return c;
        }

        public static double[,] LinearCombination(double scalar1, double[,] matrix1, double scalar2, double[,] matrix2)
        {
            int m = matrix1.GetLength(0);
            int n = matrix1.GetLength(1);
            if ((matrix2.GetLength(0) != m) || (matrix2.GetLength(1) != n)) throw new NonMatchingDimensionsException(
                "Cannot add arrays with different length");
            var c = new double[m, n];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    c[i, j] = scalar1 * matrix1[i, j] + scalar2 * matrix2[i, j];
                }
            }
            return c;
        }

        public static double DotProduct(double[] a, double[] b)
        {
            int n = a.Length;
            if (b.Length != n) throw new NonMatchingDimensionsException("Cannot add arrays with different length");
            double sum = 0.0;
            for (int i = 0; i < n; ++i)
            {
                sum += a[i] * b[i];
            }
            return sum;
        }

        public static double[] MatrixTimesVector(double[,] matrix, double[] vector)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            if (vector.Length != n) throw new NonMatchingDimensionsException("Invalid dimensions");
            double[] result = new double[m];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    result[i] += matrix[i, j] * vector[j];
                }
            }
            return result;
        }

        public static double[,] MatrixTimesMatrix(double[,] leftMatrix, double[,] rightMatrix)
        {
            if (leftMatrix.GetLength(1) != rightMatrix.GetLength(0))
                throw new NonMatchingDimensionsException("Invalid dimensions");
            double[,] result = new double[leftMatrix.GetLength(0), rightMatrix.GetLength(1)];
            for (int i = 0; i < leftMatrix.GetLength(0); ++i)
            {
                for (int j = 0; j < rightMatrix.GetLength(1); ++j)
                {
                    double sum = 0.0;
                    for (int k = 0; k < leftMatrix.GetLength(1); ++k)
                    {
                        sum += leftMatrix[i, k] * rightMatrix[k, j];
                    }
                    result[i, j] = sum;
                }
            }
            return result;
        }

        public static double[] Round(double[] vector, int decimals)
        {
            int n = vector.Length;
            double[] rounded = new double[n];
            for (int i = 0; i < n; ++i)
            {
                rounded[i] = Math.Round(vector[i], decimals);
            }
            return rounded;
        }

        public static double[,] Round(double[,] matrix, int decimals)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            double[,] rounded = new double[m, n];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    rounded[i, j] = Math.Round(matrix[i, j], decimals);
                }
            }
            return rounded;
        }

        public static double[] Scale(double scalar, double[] vector)
        {
            int n = vector.Length;
            double[] scaled = new double[n];
            for (int i = 0; i < n; ++i)
            {
                scaled[i] = scalar * vector[i];
            }
            return scaled;
        }

        public static double[,] Scale(double scalar, double[,] matrix)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            double[,] scaled = new double[m, n];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    scaled[i,j] = scalar * matrix[i, j];
                }
            }
            return scaled;
        }

        public static double[,] Transpose(double[,] matrix)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            double[,] transpose = new double[n, m];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    transpose[j, i] = matrix[i, j];
                }
            }
            return transpose;
        }
    }
}
