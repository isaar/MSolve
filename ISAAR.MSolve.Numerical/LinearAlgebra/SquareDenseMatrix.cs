using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class SquareDenseMatrix: DenseMatrix, ISquareMatrix
    {
        private static readonly double zeroTolerance = 1e-15;

        public int Order { get; }

        #region construction
        private SquareDenseMatrix(double[,] data) : base(data)
        {
            this.Order = data.GetLength(0);
        }

        public static SquareDenseMatrix CreateZero(int order)
        {
            return new SquareDenseMatrix(new double[order, order]);
        }

        public static SquareDenseMatrix CreateWithValue(int order, double value)
        {
            double[,] data = new double[order, order];
            for (int r = 0; r < order; ++r)
            {
                for (int c = 0; c < order; ++c)
                {
                    data[r, c] = value;
                }
            }
            return new SquareDenseMatrix(data);
        }

        // No you stupid thing, I do not want to hide inheritted static methods. They should not be inheritted at all.
        public static new SquareDenseMatrix CreateFromArray(double[,] data, bool copyArray = true)
        {
            if (data.GetLength(0) != data.GetLength(1))
            {
                string msg = string.Format("Provided array must have the same dimensions, but was ({0}x{1})",
                    data.GetLength(0), data.GetLength(1));
                throw new NonMatchingDimensionsException(msg);
            }
            if (copyArray)
            {
                double[,] clone = new double[data.GetLength(0), data.GetLength(1)];
                for (int r = 0; r < data.GetLength(0); ++r)
                {
                    for (int c = 0; c < data.GetLength(1); ++c)
                    {
                        clone[r, c] = data[r, c];
                    }
                }
                return new SquareDenseMatrix(clone);
            }
            else return new SquareDenseMatrix(data);
        }

        public static SquareDenseMatrix CreateFromMatrix(SquareDenseMatrix matrix)
        {
            return new SquareDenseMatrix(matrix.data);
        }
        #endregion

        public double CalcDeterminant()
        {
            try
            {
                return FactorLU(false).CalcDeterminant();
            }
            catch (SingularMatrixException ex)
            {
                return 0.0;
            }
        }

        public ILUFactorization FactorLU(bool ldu = false)
        {
            return SquareDenseMatrix.CreateFromMatrix(this).FactorLUIntoThis(ldu);
        }

        public ILUFactorization FactorLUIntoThis(bool ldu = false)
        {
            if (ldu == true) throw new NotImplementedException();

            int[] permutation = new int[Order];

            // Index k = current column to work for. The entries of U below the diagonal will be made 0.
            // Index i = row index. Only rows i >= k will be processed
            // Index j = column index. Only columns j >= k will be processed.
            for (int k = 0; k < Order; ++k) 
            {
                // Partial pivoting: find the max |A[i, k]| for i>=k
                int pivotRow = k;
                double max = Math.Abs(data[k, k]);
                for (int i = k + 1; i < Order; ++i)
                {
                    double abs = Math.Abs(data[i, k]);
                    if (abs > max)
                    {
                        max = abs;
                        pivotRow = i;
                    }
                }
                if (max <= zeroTolerance)
                {
                    string msg = string.Format("Zero pivot detected during (outer) iteration {0} of elimination", k);
                    throw new SingularMatrixException(msg);
                }

                // Partial pivoting: exhange the current row with pivotRow in L and U
                permutation[k] = pivotRow;
                permutation[pivotRow] = k;
                double swap;
                for (int j = 0; j < Order; ++j)
                {
                    swap = data[k, j];
                    data[k, j] = data[pivotRow, j];
                    data[pivotRow, j] = swap;
                }

                // Row reduction for rows i > k
                for (int i = k + 1; i < Order; ++i)
                {
                    double coeff = data[i, k] / data[k, k];
                    data[i, k] = coeff; // Store the row reduction coefficient in the lower triangle
                    for (int j = k; j < Order; ++j) // Reduce: row i  = row i - coeff * row k
                    {
                        data[i, j] -= coeff * data[k, j];
                    }
                }
            }

            return new LUFactorization(data, permutation);
        }

        /// <summary>
        /// Creates the inverse of this matrix using Gauss-Jordan. 
        /// Throws SingularMatrixException if Gauss elimination fails due to zero pivots.
        /// For an nxn matrix, the complexity is O(n^3) and the an extra 2*n^2 memory is needed. 
        /// </summary>
        /// <returns></returns>
        public ISquareMatrix Invert()
        {
            double[,] matrix = this.CopyToArray(); // We need a copy to modify during elimination
            double[,] inverse = new double[Order, Order];
            for (int i = 0; i < Order; ++i) inverse[i, i] = 1.0; //Begin with the identity matrix;

            // Index k = current column to work for. The entries of U below the diagonal will be made 0.
            // Index i = row index. Only rows i >= k will be processed
            // Index j = column index. Only columns j >= k will be processed.
            for (int k = 0; k < Order; ++k)
            {
                // Partial pivoting: find the max |A[i, k]| for i>=k
                int pivotRow = k;
                double max = Math.Abs(data[k, k]);
                for (int i = k + 1; i < Order; ++i)
                {
                    double abs = Math.Abs(data[i, k]);
                    if (abs > max)
                    {
                        max = abs;
                        pivotRow = i;
                    }
                }
                if (max <= zeroTolerance)
                {
                    string msg = string.Format("Zero pivot detected during (outer) iteration {0} of elimination", k);
                    throw new SingularMatrixException(msg);
                }

                // Partial pivoting: exhange the current row with pivotRow in the matrix and its inverse
                double swap;
                for (int j = 0; j < Order; ++j)
                {
                    // Original matrix
                    swap = matrix[k, j];
                    matrix[k, j] = matrix[pivotRow, j];
                    matrix[pivotRow, j] = swap;
                    // Inverse
                    swap = inverse[k, j];
                    inverse[k, j] = inverse[pivotRow, j];
                    inverse[pivotRow, j] = swap;
                }

                // Turn entry at pivot position to 1
                double pivot = matrix[k, k]; // It cannot be 0, since an exception would have been thrown previously 
                for (int j = 0; j < Order; ++j) // Divide row k of the original matrix with pivot. For j<k there are 0s
                {
                    matrix[k, j] /= pivot; // TODO: Optimize this: for j < k there are 0s in original matrix
                    inverse[k, j] /= pivot; // TODO: Optimize this, if possible
                }

                // Row reduction
                for (int i = 0; i < Order; ++i)
                {
                    if (i == k) continue;
                    double coeff = matrix[i, k]; // The pivot is now 1 so there is no need to divide by it.
                    for (int j = 0; j < Order; ++j) // Reduce: row i  = row i - coeff * row k
                    {
                        matrix[i, j] -= coeff * matrix[k, j];
                        inverse[i, j] -= coeff * inverse[k, j];
                        // TODO: Optimizations for original matrix:
                        // 1) For i>k, only j>=k need to be accessed. 
                        // 2) After all i>k are completed, then i<k can be performed in reverse order.
                        // 3) Then for each diagonal k, row k is the row of I. Thus reducing the column above it
                        // can be done by just setting them to 0, since 
                        // coeff = matrix[i,k] / matrix[k,k] = matrix[i,k] and for j!=k : 
                        // matrix[i,j] = matrix[i,j] - matrix[k,j]*coeff = matrix[i,j] - 0*coeff = matrix[k,j]
                        // , while for j=k:
                        // matrix[i,k] = matrix[i,k] - matrix[k,k]*coeff = matrix[i,k] - 1*matrix[i,k] = 0
                        // 4) However coeff must be used to reduce the rows of the inverse matrix as usual. Not sure if 
                        // there is a certain pattern of zeros there.
                    }
                }
            }

            return new SquareDenseMatrix(inverse); // The copy of the original matrix is destroyed after this point.
        }
    }
}
