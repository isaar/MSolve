using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities
{
    public static class MatrixOperations
    {
        public static double[] Add(double[] a, double[] b)
        {
            int n = a.Length;
            if (b.Length != n) throw new RankException("Cannot add arrays with different length");
            var c = new double[n];
            for (int i = 0; i < n; ++i)
            {
                c[i] = a[i] + b[i];
            }
            return c;
        }

        public static double[,] Add(double[,] a, double[,] b)
        {
            int m = a.GetLength(0);
            int n = a.GetLength(1);
            if ((b.GetLength(0) != m) || (b.GetLength(1) != n)) throw new RankException(
                "Cannot add arrays with different length");
            var c = new double[m, n];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < 0; ++j)
                {
                    c[i, j] = a[i, j] + b[i, j];
                }
            }
            return c;
        }

        public static double DotProduct(double[] a, double[] b)
        {
            int n = a.Length;
            if (b.Length != n) throw new RankException("Cannot add arrays with different length");
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
            if (vector.Length != n) throw new ArgumentException("Invalid dimensions");
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
            for (int i = 0; i < n; ++i)
            {
                vector[i] *= scalar;
            }
            return vector;
        }

        public static double[,] Scale(double scalar, double[,] matrix)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < 0; ++j)
                {
                    matrix[i, j] *= scalar;
                }
            }
            return matrix;
        }
    }
}
