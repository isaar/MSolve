using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: this should not be in the same folder with the actual matrices and their interfaces
//TODO: Split this into many classes: IMatrixExtensions, MatrixExtensions, CsrMatrixExtensions (or collectively 
//      SparseMatrixExtensions, etc).
//TODO: The GetColumn, GetDiagonal, GetRow should be implemented as default interface methods and overwritten by concrete matrix  
//      classes when possible.
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Defines common matrix operations that can be used as extensions for methods for matrix classes or interfaces. 
    /// If one of the following methods is also implemented in a concrete class, use that one, as it will be far more efficient.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MatrixExtensions
    {
        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="matrix1"/>[i, j] + <paramref name="matrix2"/>[i, j], 
        /// for 0 &lt;= i &lt; <see cref="IIndexable2D.NumRows"/>, 0 &lt;= j &lt; <see cref="IIndexable2D.NumColumns"/>.
        /// The resulting entries are written to a new <see cref="IMatrixView"/> instance.
        /// </summary>
        /// <param name="matrix1">The first <see cref="IMatrixView"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix2"/>.</param>
        /// <param name="matrix2">The second <see cref="IMatrixView"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix1"/> and <paramref name="matrix2"/>
        ///     have a different number of <see cref="IIndexable2D.NumRows"/> or 
        ///     <see cref="IIndexable2D.NumColumns"/>.</exception>
        public static IMatrix Add(this IMatrixView matrix1, IMatrixView matrix2) => matrix1.Axpy(matrix2, 1.0);

        /// <summary>
        /// Performs the operation: 
        /// <paramref name="matrix1"/>[i, j] = <paramref name="matrix1"/>[i, j] + <paramref name="matrix2"/>[i, j], 
        /// for 0 &lt;= i &lt; <see cref="IIndexable2D.NumRows"/>, 0 &lt;= j &lt; <see cref="IIndexable2D.NumColumns"/>.
        /// The resulting matrix overwrites the entries of <paramref name="matrix1"/>.
        /// </summary>
        /// <param name="matrix1">The first <see cref="IMatrix"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix2"/>.</param>
        /// <param name="matrix2">The second <see cref="IMatrixView"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix1"/> and <paramref name="matrix2"/>
        ///     have a different number of <see cref="IIndexable2D.NumRows"/> or 
        ///     <see cref="IIndexable2D.NumColumns"/>.</exception>
        /// <exception cref="PatternModifiedException">Thrown if an <paramref name="matrix1"/>[i, j] needs to be 
        ///     overwritten, but that is not permitted by the matrix storage format.</exception>
        public static void AddIntoThis(this IMatrix matrix1, IMatrixView matrix2) => matrix1.AxpyIntoThis(matrix2, 1.0);

        /// <summary>
        /// Iterates over the non-zero entries of the matrix, which are defined as the entries such that: 
        /// <paramref name="matrix"/>[i,j].
        /// </summary>
        /// <param name="matrix">The matrix whose non-zero entries will be iterated over.</param>
        public static IEnumerable<(int row, int col, double val)> EnumerateNonZeros(this IIndexable2D matrix)
        {
            for (int j = 0; j < matrix.NumColumns; ++j)
            {
                for (int i = 0; i < matrix.NumRows; ++i)
                {
                    if (matrix[i, j] != 0.0) yield return (i, j, matrix[i, j]);
                }
            }
        }

        /// <summary>
        /// Iterates over the non-zero entries of the matrix, which are defined as the entries such that: 
        /// Math.Abs(<paramref name="matrix"/>[i,j]) &gt; 0.
        /// </summary>
        /// <param name="matrix">The matrix whose non-zero entries will be iterated over.</param>
        /// <param name="tolerance">A small value under which an entry is considered to be 0.</param>
        public static IEnumerable<(int row, int col, double val)> EnumerateNonZeros(this IIndexable2D matrix, double tolerance)
        {
            for (int j = 0; j < matrix.NumColumns; ++j)
            {
                for (int i = 0; i < matrix.NumRows; ++i)
                {
                    if (Math.Abs(matrix[i, j]) > 0.0) yield return (i, j, matrix[i, j]);
                }
            }
        }

        /// <summary>
        /// Returns a <see cref="Vector"/> with the entries of the matrix's main diagonal.
        /// </summary>
        /// <param name="matrix">The matrix whose diagonal will be copied. It must be square.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        public static Vector GetDiagonal(this IMatrixView matrix) => Vector.CreateFromArray(matrix.GetDiagonalAsArray(), false);

        /// <summary>
        /// Returns an array with the entries of the matrix's main diagonal.
        /// </summary>
        /// <param name="matrix">The matrix whose diagonal will be copied. It must be square.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        public static double[] GetDiagonalAsArray(this IMatrixView matrix)
        {
            Preconditions.CheckSquare(matrix);
            double[] diag = new double[matrix.NumRows];
            for (int i = 0; i < matrix.NumRows; ++i) diag[i] = matrix[i, i];
            return diag;
        }

        /// <summary>
        /// Returns true if <paramref name="matrix"/>[i, j] and <paramref name="matrix"/>[j, i] are equal or at least within the 
        /// specified <paramref name="tolerance"/> for all 0 &lt;= i &lt; <see cref="IIndexable2D.NumRows"/>, 
        /// 0 &lt;= j &lt; <see cref="IIndexable2D.NumColumns"/>. 
        /// </summary>
        /// <param name="matrix">The matrix that will be checked for symmetry.</param>
        /// <param name="tolerance">The entries at (i, j), (j, i) the matrix will be considered equal, if
        ///     (<paramref name="matrix"/>[i, j] - <paramref name="matrix"/>[i, j]) / <paramref name="matrix"/>[i, j] 
        ///         &lt;= <paramref name="tolerance"/>. 
        ///     Setting <paramref name="tolerance"/> = 0, will check if these entries are exactly the same.</param>
        public static bool IsSymmetric(this IIndexable2D matrix, double tolerance = double.Epsilon)
        {
            var comparer = new ValueComparer(tolerance);
            if (matrix.NumRows != matrix.NumColumns) return false;
            for (int i = 0; i < matrix.NumRows; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    if (!comparer.AreEqual(matrix[i, j], matrix[j, i])) return false;
                }
            }
            return true;
        }
        
        /// <summary>
        /// Computes the Reduced Row Echelon Form (rref) of the matrix and finds the independent columns of the matrix.  
        /// See https://en.wikipedia.org/wiki/Row_echelon_form#Reduced_row_echelon_form.
        /// </param>
        public static (Matrix rref, List<int> independentCols) ReducedRowEchelonForm(this IMatrixView matrix)
            => ReducedRowEchelonForm(matrix,
                GlobalConstants.MachinePrecisionDouble * Math.Max(matrix.NumRows, matrix.NumColumns) * matrix.NormInf());

        /// <summary>
        /// Computes the Reduced Row Echelon Form (rref) of the matrix and finds the independent columns of the matrix.
        /// See https://en.wikipedia.org/wiki/Row_echelon_form#Reduced_row_echelon_form. 
        /// </summary>
        /// <param name="pivotTolerance">
        /// If the absolute values of a diagonal entry is less than this tolerance it is assumed to be zero.
        /// </param>
        public static (Matrix rref, List<int> independentCols) ReducedRowEchelonForm(this IMatrixView matrix, 
            double pivotTolerance)
        {
            // Ported from octave's built-in implementation: https://searchcode.com/codesearch/view/9591940/.
            var rref = matrix.CopyToFullMatrix();
            int numRows = rref.NumRows;
            int numCols = rref.NumColumns;

            var independentCols = new List<int>();
            int row = 0;
            for (int col = 0; col < numCols; ++col)
            {
                // Find the pivot row
                int pivotIdx = int.MinValue;
                double pivotValue = double.MinValue;
                for (int i = row; i < numRows; ++i)
                {
                    double abs = Math.Abs(rref[i, col]);
                    if (abs > pivotValue)
                    {
                        pivotIdx = i;
                        pivotValue = abs;
                    }
                }

                if (pivotValue <= pivotTolerance) 
                {
                    // Skip column c, making sure the approximately zero terms are actually zero.
                    for (int i = row; i < numRows; ++i) rref[i, col] = 0.0;
                }
                else
                {
                    // Keep track of bound variables
                    independentCols.Add(col);

                    // Swap current row and pivot row if necessary
                    if (pivotIdx != row)
                    {
                        for (int j = col; j < numCols; ++j)
                        {
                            double swap = rref[row, j];
                            rref[row, j] = rref[pivotIdx, j];
                            rref[pivotIdx, j] = swap;
                        }
                    }

                    // Normalize pivot row
                    double diagonal = rref[row, col];
                    for (int j = col; j < numCols; ++j) rref[row, j] /= diagonal;

                    // Eliminate the current column
                    for (int i = 0; i < numRows; ++i)
                    {
                        if (i == row) continue;
                        double scale = rref[i, col];
                        for (int j = col; j < numCols; ++j) rref[i, j] -= scale * rref[row, j];
                    }

                    // Check if done
                    if (row == numRows - 1) return (rref, independentCols);
                    else ++row;
                }
            }

            return (rref, independentCols);
        }

        /// <summary>
        /// Creates a new <see cref="Matrix"/> that contains the entries of <paramref name="matrix"/> with a different order,
        /// which is specified by the provided <paramref name="permutation"/> and <paramref name="oldToNew"/>.
        /// </summary>
        /// <param name="matrix">The matrix whose rows and columns will be reordered.</param>
        /// <param name="permutation">An array that contains the row/column indices of <paramref name="matrix"/> in a 
        ///     different order.</param>
        /// <param name="oldToNew">If true, 
        ///     reordered[<paramref name="permutation"/>[i], <paramref name="permutation"/>[j]] =  original[i, j]. If false, 
        ///     reordered[i, j] = original[<paramref name="permutation"/>[i], <paramref name="permutation"/>[j]].</param>
        public static Matrix Reorder(this IIndexable2D matrix, IReadOnlyList<int> permutation, bool oldToNew)
        {
            int order = matrix.NumRows;
            if (matrix.NumColumns != order)
            {
                throw new NonMatchingDimensionsException("This operation works on square matrices only.");
            }
            if (permutation.Count != order)
            {
                throw new NonMatchingDimensionsException($"This matrix has order = {order}, while the permutation vector"
                    + $" has {permutation.Count} entries.");
            }
            var reordered = Matrix.CreateZero(order, order);
            if (oldToNew)
            {
                for (int j = 0; j < order; ++j)
                {
                    int newCol = permutation[j];
                    for (int i = 0; i < order; ++i) reordered[permutation[i], newCol] = matrix[i, j];
                }
            }
            else
            {
                for (int j = 0; j < order; ++j)
                {
                    int oldCol = permutation[j];
                    for (int i = 0; i < order; ++i) reordered[i, j] = matrix[permutation[i], oldCol];
                }
            }
            return reordered;
        }

        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="matrix1"/>[i, j] - <paramref name="matrix2"/>[i, j], 
        /// for 0 &lt;= i &lt; <see cref="IIndexable2D.NumRows"/>, 0 &lt;= j &lt; <see cref="IIndexable2D.NumColumns"/>.
        /// The resulting entries are written to a new <see cref="IMatrixView"/> instance.
        /// </summary>
        /// <param name="matrix1">The first <see cref="IMatrixView"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix2"/>.</param>
        /// <param name="matrix2">The second <see cref="IMatrixView"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix1"/> and <paramref name="matrix2"/>
        ///     have a different number of <see cref="IIndexable2D.NumRows"/> or 
        ///     <see cref="IIndexable2D.NumColumns"/>.</exception>
        public static IMatrix Subtract(this IMatrixView matrix1, IMatrixView matrix2) => matrix1.Axpy(matrix2, -1.0);

        /// <summary>
        /// Performs the operation: 
        /// <paramref name="matrix1"/>[i, j] = <paramref name="matrix1"/>[i, j] + <paramref name="matrix2"/>[i, j], 
        /// for 0 &lt;= i &lt; <see cref="IIndexable2D.NumRows"/>, 0 &lt;= j &lt; <see cref="IIndexable2D.NumColumns"/>.
        /// The resulting matrix overwrites the entries of <paramref name="matrix1"/>.
        /// </summary>
        /// <param name="matrix1">The first <see cref="IMatrix"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix2"/>.</param>
        /// <param name="matrix2">The second <see cref="IMatrixView"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix1"/> and <paramref name="matrix2"/>
        ///     have a different number of <see cref="IIndexable2D.NumRows"/> or 
        ///     <see cref="IIndexable2D.NumColumns"/>.</exception>
        /// <exception cref="PatternModifiedException">Thrown if an <paramref name="matrix1"/>[i, j] needs to be 
        ///     overwritten, but that is not permitted by the matrix storage format.</exception>
        public static void SubtractIntoThis(this IMatrix matrix1, IMatrixView matrix2) => matrix1.AxpyIntoThis(matrix2, -1.0);
    }
}
