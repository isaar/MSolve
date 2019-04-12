using System;
using CSparse;
using CSparse.Double;
using CSparse.Double.Factorization;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Implement IIndexable2D to allow easy output.
//TODO: Allow other orderings, as I do in CholeskySuiteSparse
//TODO: CSparse.NET also provides matrix update and downdate operations.
//TODO: Improve error checking
//TODO: expose internal factorized arrays 
namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// LU factorization of a sparse square matrix using the CSparse.NET library. The original matrix must be in 
    /// Compressed Sparse Columns format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class LUCSparseNet : ITriangulation
    {
        public const double defaultPivotTolelance = 0.1; // between 0.0, 1.0. TODO: Find a good default
        private readonly SparseLU factorization;

        private LUCSparseNet(int order, SparseLU factorization)
        {
            this.Order = order;
            this.factorization = factorization;
        }

        /// <summary>
        /// The number of non-zero entries (and explicitly stored zeros) in the explicitly stored lower and upper triangular 
        /// factors after LU factorization.
        /// </summary>
        public int NumNonZerosUpper => factorization.NonZerosCount;

        /// <summary>
        /// The number of rows/columns of the square matrix. 
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Performs the LU factorization: A = L * U of a ssquare matrix A.  The matrix A is provided in CSC format by 
        /// <paramref name="cscValues"/>, <paramref name="cscRowIndices"/> and <paramref name="cscColOffsets"/>. 
        /// </summary>
        /// <param name="order">The number of rows/columns of the square matrix.</param>
        /// <param name="numNonZeros">The number of explicitly stored entries of the matrix before the factorization.</param>
        /// <param name="cscValues">
        /// Contains the non-zero entries of the upper triangle. Its length must be equal to <paramref name="numNonZeros"/>.
        /// The non-zero entries of each row must appear consecutively in <paramref name="cscValues"/>. They should also be 
        /// sorted in increasing order of their row indices, to speed up subsequent the factorization. 
        /// </param>
        /// <param name="cscRowIndices">
        /// Contains the row indices of the non-zero entries. Its length must be equal to <paramref name="numNonZeros"/>. 
        /// There is an 1 to 1 matching between these two arrays: <paramref name="cscRowIndices"/>[i] is the row index of the 
        /// entry <paramref name="cscValues"/>[i]. Also: 0 &lt;= <paramref name="cscRowIndices"/>[i] &lt; 
        /// <paramref name="order"/>.
        /// </param>
        /// <param name="cscColOffsets">
        /// Contains the index of the first entry of each column into the arrays <paramref name="cscValues"/> and 
        /// <paramref name="cscRowIndices"/>. Its length must be <paramref name="order"/> + 1. The last entry must be 
        /// <paramref name="numNonZeros"/>.
        /// </param>
        /// <param name="pivotTolerance">The partial pivoting tolerance (from 0.0 to 1.0).</param>
        /// <exception cref="SingularMatrixException">Thrown if the original matrix is not positive definite.</exception>
        public static LUCSparseNet Factorize(int order, int numNonZeros, double[] cscValues, int[] cscRowIndices,
            int[] cscColOffsets, double pivotTolerance = defaultPivotTolelance)
        {
            try
            {
                var matrixCSparse = new SparseMatrix(order, order, cscValues, cscRowIndices, cscColOffsets);
                var factorization = SparseLU.Create(matrixCSparse, ColumnOrdering.Natural, pivotTolerance);
                return new LUCSparseNet(order, factorization);
            }
            catch (Exception ex) //TODO: how can I make sure this exception was thrown because of an indefinite matrix?
            {
                throw new SingularMatrixException(ex.Message);
            }
        }

        /// <summary>
        /// Performs the LU factorization: A = L * U of a square matrix A. The matrix A is provided in CSC format. 
        /// </summary>
        /// <param name="matrix">The matrix in symmetric (only upper triangle) CSC format.</param>
        /// <param name="pivotTolerance">The partial pivoting tolerance (from 0.0 to 1.0).</param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the original matrix is not positive definite.</exception>
        public static LUCSparseNet Factorize(CscMatrix matrix, double pivotTolerance = defaultPivotTolelance)
            => Factorize(matrix.NumColumns, matrix.NumNonZeros, matrix.RawValues, matrix.RawRowIndices, 
                matrix.RawColOffsets, pivotTolerance);

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        public double CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// See <see cref="ITriangulation.SolveLinearSystem(Vector, Vector)"/>.
        /// </summary>
        public void SolveLinearSystem(Vectors.Vector rhs, Vectors.Vector solution)
            => factorization.Solve(rhs.RawData, solution.RawData);
    }
}
