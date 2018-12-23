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
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    /// <summary>
    /// Cholesky factorization of a sparse symmetric positive definite matrix using the CSparse.NET library. The original matrix 
    /// must be in Compressed Sparse Columns format, with only the upper triangle stored. This class may serve as a managed
    /// alternative to <see cref="CholeskySuiteSparse"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CholeskyCSparseNet : ITriangulation
    {
        private readonly double[] cscValues;
        private readonly int[] cscRowIndices, cscColOffsets;
        private readonly SparseCholesky factorization;

        private CholeskyCSparseNet(int order, double[] cscValues, int[] cscRowIndices,
            int[] cscColOffsets, SparseCholesky factorization)
        {
            this.Order = order;
            this.cscValues = cscValues;
            this.cscRowIndices = cscRowIndices;
            this.cscColOffsets = cscColOffsets;
            this.factorization = factorization;
        }

        /// <summary>
        /// The number of non-zero entries (and explicitly stored zeros) in the explicitly stored upper triangular factor 
        /// after Cholesky factorization.
        /// </summary>
        public int NumNonZerosUpper => factorization.NonZerosCount;

        /// <summary>
        /// The number of rows/columns of the square matrix. 
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Performs the Cholesky factorization: A = L * L^T of a symmetric positive definite matrix A. 
        /// Only the upper triangle of the original matrix is required and is provided in symmetric CSC format by 
        /// <paramref name="cscValues"/>, <paramref name="cscRowIndices"/> and <paramref name="cscColOffsets"/>. 
        /// The user may choose between supernodal or simplicial factorization. It is also possible to automatically reorder 
        /// the matrix, using the algorithms provided by SuiteSparse.
        /// The factorized data, which may be sufficiently larger than the original matrix due to fill-in, will be written to 
        /// unmanaged memory.
        /// </summary>
        /// <param name="order">The number of rows/columns of the square matrix.</param>
        /// <param name="numNonZerosUpper">The number of explicitly stored entries in the upper triangle of the matrix.</param>
        /// <param name="cscValues">
        /// Contains the non-zero entries of the upper triangle. Its length must be equal to <paramref name="numNonZerosUpper"/>.
        /// The non-zero entries of each row must appear consecutively in <paramref name="cscValues"/>. They should also be 
        /// sorted in increasing order of their row indices, to speed up subsequent the factorization. 
        /// </param>
        /// <param name="cscRowIndices">
        /// Contains the row indices of the non-zero entries. Its length must be equal to <paramref name="numNonZerosUpper"/>. 
        /// There is an 1 to 1 matching between these two arrays: <paramref name="cscRowIndices"/>[i] is the row index of the 
        /// entry <paramref name="cscValues"/>[i]. Also: 0 &lt;= <paramref name="cscRowIndices"/>[i] &lt; 
        /// <paramref name="order"/>.
        /// </param>
        /// <param name="cscColOffsets">
        /// Contains the index of the first entry of each column into the arrays <paramref name="cscValues"/> and 
        /// <paramref name="cscRowIndices"/>. Its length must be <paramref name="order"/> + 1. The last entry must be 
        /// <paramref name="numNonZerosUpper"/>.
        /// </param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the original matrix is not positive definite.</exception>
        public static CholeskyCSparseNet Factorize(int order, int numNonZerosUpper, double[] cscValues, int[] cscRowIndices,
            int[] cscColOffsets)
        {
            try
            {
                var matrixCSparse = new SparseMatrix(order, order, cscValues, cscRowIndices, cscColOffsets);
                var factorization = SparseCholesky.Create(matrixCSparse, ColumnOrdering.Natural);
                return new CholeskyCSparseNet(order, cscValues, cscRowIndices, cscColOffsets, factorization);
            }
            catch (Exception ex) //TODO: how can I make sure this exception was thrown because of an indefinite matrix?
            {
                throw new IndefiniteMatrixException(ex.Message);
            }
        }

        /// <summary>
        /// Performs the Cholesky factorization: A = L * L^T of a symmetric positive definite matrix A. 
        /// Only the upper triangle of the original matrix is required and is provided in symmetric CSC format. 
        /// The user may choose between supernodal or simplicial factorization. It is also possible to automatically reorder 
        /// the matrix, using the algorithms provided by SuiteSparse.
        /// The factorized data, which may be sufficiently larger than the original matrix due to fill-in, will be written to 
        /// unmanaged memory.
        /// </summary>
        /// <param name="matrix">The matrix in symmetric (only upper triangle) CSC format.</param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the original matrix is not positive definite.</exception>
        public static CholeskyCSparseNet Factorize(SymmetricCscMatrix matrix)
            => Factorize(matrix.NumColumns, matrix.NumNonZerosUpper, matrix.RawValues, matrix.RawRowIndices, 
                matrix.RawColOffsets);

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
