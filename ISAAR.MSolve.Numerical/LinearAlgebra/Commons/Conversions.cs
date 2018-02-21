using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Commons
{
    public static class Conversions
    {
        public static double[,] Array2DLowerToSymmetric(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[,] symm = new double[n, n];
            Array.Copy(array2D, symm, n * n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    symm[j, i] = symm[i, j];
                }
            }
            return symm;
        }

        public static double[,] Array2DUpperToSymmetric(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[,] symm = new double[n, n];
            Array.Copy(array2D, symm, n * n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    symm[i, j] = symm[j, i];
                }
            }
            return symm;
        }

        public static double[] Array2DToFullColMajor(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numColumns = array2D.GetLength(1);
            double[] array1D = new double[numRows * numColumns];
            for (int j = 0; j < numColumns; ++j)
            {
                for (int i = 0; i < numRows; ++i)
                {
                    array1D[j * numRows + i] = array2D[i, j];
                }
            }
            return array1D;
        }

        public static double[] Array2DToFullRowMajor(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numColumns = array2D.GetLength(1);
            double[] array1D = new double[numRows * numColumns];
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numColumns; ++j)
                {
                    array1D[i * numColumns + j] = array2D[i, j];
                }
            }
            return array1D;
        }

        public static double[] Array2DToPackedLowerColMajor(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[] array1D = new double[(n * (n + 1)) / 2];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int j = 0; j < n; ++j)
            {
                for (int i = j; i < n; ++i)
                {
                    array1D[counter] = array2D[i, j];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array1D;
        }

        public static double[] Array2DToPackedLowerRowMajor(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[] array1D = new double[(n * (n + 1)) / 2];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    array1D[counter] = array2D[i, j];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array1D;
        }

        public static double[] Array2DToPackedUpperColMajor(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[] array1D = new double[(n * (n + 1)) / 2];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    array1D[counter] = array2D[i, j];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array1D;
        }

        public static double[] Array2DToPackedUpperRowMajor(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[] array1D = new double[(n * (n + 1)) / 2];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int i = 0; i < n; ++i)
            {
                for (int j = i; j < n; ++j)
                {
                    array1D[counter] = array2D[i, j];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array1D;
        }

        public static double[,] FullColMajorToArray2D(double[] array1D, int numRows, int numColumns)
        {
            double[,] array2D = new double[numRows, numColumns];
            for (int j = 0; j < numColumns; ++j)
            {
                for (int i = 0; i < numRows; ++i)
                {
                    array2D[i, j] = array1D[j * numRows + i];
                }
            }
            return array2D;
        }

        /// <summary>
        /// Extract the lower part only. Both I/O are column major with full storage.
        /// </summary>
        /// <param name="full"></param>
        /// <param name="unitDiagonal"></param>
        /// <returns></returns>
        public static double[] FullColMajorToFullLowerColMajor(double[] full, bool unitDiagonal)
        {
            int n = FullLengthToOrder(full.Length);
            double[] lower = new double[n * n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = j + 1; i < n; ++i)
                {
                    lower[j * n + i] = full[j * n + i];
                }
            }
            if (unitDiagonal)
            {
                for (int d = 0; d < n; ++d) lower[d * n + d] = 1.0;
            }
            else
            {
                for (int d = 0; d < n; ++d) lower[d * n + d] = full[d * n + d];
            }
            return lower;
        }

        /// <summary>
        /// Extract the upper part only. Both I/O are column major with full storage.
        /// </summary>
        /// <param name="full"></param>
        /// <param name="unitDiagonal"></param>
        /// <returns></returns>
        public static double[] FullColMajorToFullUpperColMajor(double[] full, bool unitDiagonal)
        {
            int n = FullLengthToOrder(full.Length);
            double[] upper = new double[n * n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < j; ++i)
                {
                    upper[j * n + i] = full[j * n + i];
                }
            }
            if (unitDiagonal)
            {
                for (int d = 0; d < n; ++d) upper[d * n + d] = 1.0;
            }
            else
            {
                for (int d = 0; d < n; ++d) upper[d * n + d] = full[d * n + d];
            }
            return upper;
        }

        public static double[,] FullRowMajorToArray2D(double[] array1D, int numRows, int numColumns)
        {
            double[,] array2D = new double[numRows, numColumns];
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numColumns; ++j)
                {
                    array2D[i, j] = array1D[i * numColumns + j];
                }
            }
            return array2D;
        }

        public static double[,] FullLowerColMajorToArray2D(double[] array1D, bool unitDiagonal)
        {
            int n = FullLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = j + 1; i < n; ++i)
                {
                    array2D[i, j] = array1D[j * n + i];
                }
            }
            if (unitDiagonal)
            {
                for (int d = 0; d < n; ++d) array2D[d, d] = 1.0;
            }
            else
            {
                for (int d = 0; d < n; ++d) array2D[d, d] = array1D[d * n + d];
            }
            return array2D;
        }

        public static double[,] FullUpperColMajorToArray2D(double[] array1D, bool unitDiagonal)
        {
            int n = FullLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < j; ++i)
                {
                    array2D[i, j] = array1D[j * n + i];
                }
            }
            if (unitDiagonal)
            {
                for (int d = 0; d < n; ++d) array2D[d, d] = 1.0;
            }
            else
            {
                for (int d = 0; d < n; ++d) array2D[d, d] = array1D[d * n + d];
            }
            return array2D;
        }

        public static double[,] PackedLowerColMajorToArray2D(double[] array1D)
        {
            int n = PackedLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int j = 0; j < n; ++j)
            {
                for (int i = j; i < n; ++i)
                {
                    array2D[i, j] = array1D[counter];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array2D;
        }

        public static double[,] PackedLowerRowMajorToArray2D(double[] array1D)
        {
            int n = PackedLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    array2D[i, j] = array1D[counter];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array2D;
        }

        public static double[,] PackedUpperColMajorToArray2D(double[] array1D)
        {
            int n = PackedLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    array2D[i, j] = array1D[counter];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array2D;
        }

        public static double[,] PackedUpperRowMajorToArray2D(double[] array1D)
        {
            int n = PackedLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int i = 0; i < n; ++i)
            {
                for (int j = i; j < n; ++j)
                {
                    array2D[i, j] = array1D[counter];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array2D;
        }

        public static int FullLengthToOrder(int length)
        {
            // length = n^2 => n = sqrt(length)
            double n = Math.Sqrt(length);
            int order = (int)Math.Round(n);
            if (order * order != length)
            {
                throw new ArgumentException("The length of the 1D array must be an integer L such that L=n*n");
            }
            return order;
        }

        public static int PackedLengthToOrder(int length)
        {
            // length = n*(n+1)/2 => n = ( -1+sqrt(1+8*length) )/2
            double n = (-1.0 + Math.Sqrt(1 + 8 * length)) / 2;
            int order = (int)Math.Round(n);
            if (order * (order + 1) / 2 != length)
            {
                throw new ArgumentException("The length of the 1D array must be an integer L such that L=n*(n+1)/2");
            }
            return order;
        }
    }
}
