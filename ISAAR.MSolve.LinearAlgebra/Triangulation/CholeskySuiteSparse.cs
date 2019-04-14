using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Providers.PInvoke;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: SuiteSparse Common should be represented here by an IDisposable class SuiteSparseCommon.
//TODO: Perhaps I should express the back/forward/full solve using the L*D*L^T, L*L^T, L^T, L^T*D factors as in CHOLMOD.
//TODO: During a non linear or dynamic analysis the sparsity pattern stays the same. Perhaps the symbolic factorization phase 
//      (and anything else that depends only on the pattern) can be cached.
//TODO: the user should be able to choose a reordering algorithm or try them all, independently of the factorization. He could 
//      then use that reordering during other system solutions with the same dofs (e.g. during a non linear analysis). As it is
//      now, the option to try various reorderings is only possible during factorization. Thus, during each linear system 
//      solution, the rhs and solution vectors must be permuted and the same pattern must be processed to find a good reordering.
namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// Cholesky factorization of a sparse symmetric positive definite matrix using the SuiteSparse library. The original matrix 
    /// must be in Compressed Sparse Columns format, with only the upper triangle stored. SuiteSparse is very efficient for  
    /// sparse matrices and provides a lot of functionality, but requires handling unmanaged memory, which is abstracted in this
    /// class. If SuiteSparse dlls are not available, try using the managed alternative <see cref="CholeskyCSparseNet"/> instead.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CholeskySuiteSparse : ITriangulation, IDisposable
    {
        private IntPtr common;
        private IntPtr factorizedMatrix;

        private CholeskySuiteSparse(int order, IntPtr suiteSparseCommon, IntPtr factorizedMatrix)
        {
            this.Order = order;
            this.common = suiteSparseCommon;
            this.factorizedMatrix = factorizedMatrix;
        }

        ~CholeskySuiteSparse()
        {
            ReleaseResources();
        }

        /// <summary>
        /// The number of non-zero entries (and explicitly stored zeros) in the explicitly stored upper triangular factor 
        /// after Cholesky factorization.
        /// </summary>
        public int NumNonZerosUpper { get => SuiteSparsePInvokes.GetFactorNonZeros(factorizedMatrix); }

        /// <summary>
        /// The number of rows/columns of the square matrix. 
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Performs the Cholesky factorization: A = L * L^T or A = L * D * L^T of a symmetric positive definite matrix A. 
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
        /// <param name="superNodal">
        /// If true, a supernodal factorization will be performed, which results in faster back/forward substitutions during the 
        /// linear system solution. If false, a simplicial factorization will be performed, which allows manipulating the 
        /// factorized matrix (e.g. using <see cref="AddRow(int, SparseVector)"/> or <see cref="DeleteRow(int)"/>.
        /// </param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the original matrix is not positive definite.</exception>
        /// <exception cref="SuiteSparseException">
        /// Thrown if the calls to SuiteSparse library fail. This usually happens if the SuiteSparse .dlls are not available or 
        /// if there is not sufficient memory to perform the factorization.
        /// </exception>
        public static CholeskySuiteSparse Factorize(int order, int numNonZerosUpper, double[] cscValues, int[] cscRowIndices,
            int[] cscColOffsets, bool superNodal)
        {
            int factorizationType = superNodal ? 1 : 0;
            IntPtr common = SuiteSparsePInvokes.CreateCommon(factorizationType, (int)SuiteSparseOrdering.Natural);
            if (common == IntPtr.Zero) throw new SuiteSparseException("Failed to initialize SuiteSparse.");
            int status = SuiteSparsePInvokes.FactorizeCSCUpper(order, numNonZerosUpper, cscValues, cscRowIndices, cscColOffsets,
                out IntPtr factorizedMatrix, common);
            if (status == -2)
            {
                SuiteSparsePInvokes.DestroyCommon(ref common);
                throw new SuiteSparseException("Factorization did not succeed. This could be caused by insufficent memory,"
                    + " due to excessive fill-in.");
            }
            else if (status >= 0)
            {
                SuiteSparsePInvokes.DestroyCommon(ref common);
                throw new IndefiniteMatrixException("The matrix is not positive definite."
                    + $" Cholesky failed at column {status} (0-based indexing).");
            }
            else return new CholeskySuiteSparse(order, common, factorizedMatrix);
        }

        /// <summary>
        /// Performs the Cholesky factorization: A = L * L^T or A = L * D * L^T of a symmetric positive definite matrix A. 
        /// Only the upper triangle of the original matrix is required and is provided in symmetric CSC format. 
        /// The user may choose between supernodal or simplicial factorization. It is also possible to automatically reorder 
        /// the matrix, using the algorithms provided by SuiteSparse.
        /// The factorized data, which may be sufficiently larger than the original matrix due to fill-in, will be written to 
        /// unmanaged memory.
        /// </summary>
        /// <param name="matrix">The matrix in symmetric (only upper triangle) CSC format.</param>
        /// <param name="superNodal">
        /// If true, a supernodal factorization will be performed, which results in faster back/forward substitutions during the 
        /// linear system solution. If false, a simplicial factorization will be performed, which allows manipulating the 
        /// factorized matrix (e.g. using <see cref="AddRow(int, SparseVector)"/> or <see cref="DeleteRow(int)"/>.
        /// </param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the original matrix is not positive definite.</exception>
        /// <exception cref="SuiteSparseException">
        /// Thrown if the calls to SuiteSparse library fail. This usually happens if the SuiteSparse .dlls are not available or 
        /// if there is not sufficient memory to perform the factorization.
        /// </exception>
        public static CholeskySuiteSparse Factorize(SymmetricCscMatrix matrix, bool superNodal)
            => Factorize(matrix.NumColumns, matrix.NumNonZerosUpper, matrix.RawValues, matrix.RawRowIndices, 
                matrix.RawColOffsets, false);

        /// <summary>
        /// Sets the <paramref name="rowIdx"/>-th row/column of the factorized matrix to the one it would have if 
        /// <paramref name="newRow"/> was set as the <paramref name="rowIdx"/>-th row/column of the original matrix and then the 
        /// factorization was performed. The existing <paramref name="rowIdx"/>-th row/column of the original matrix must be 
        /// equal to the <paramref name="rowIdx"/>-th row/col of the identity matrix. 
        /// </summary>
        /// <param name="rowIdx">The index of the row/column to modify. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; this.<see cref="Order"/>.</param>
        /// <param name="newRow">The entries of the row/column before factorization. Constraints:
        ///     <paramref name="newRow"/>.<see cref="IIndexable1D.Length"/> == this.<see cref="Order"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowIdx"/> violates the described constraints.
        ///     </exception>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="newRow"/> violates the described 
        ///     constraints. </exception>
        /// <exception cref="SuiteSparseException">Thrown if the call to SuiteSparse library fails, usually due to insufficient 
        ///     memory.</exception>
        public void AddRow(int rowIdx, SparseVector newRow)
        {
            //TODO: use Preconditions for these tests and implement IIndexable2D.
            if ((rowIdx < 0) || (rowIdx >= Order))
            {
                throw new IndexOutOfRangeException($"Cannot access row {rowIdx} in a"
                    + $" {Order}-by-{Order} matrix");
            }
            if (newRow.Length != Order)
            {
                throw new NonMatchingDimensionsException($"The new row/column must have the same number of rows as this"
                    + $"{Order}-by-{Order} factorized matrix, but was {newRow.Length}-by-1");
            }

            int nnz = newRow.CountNonZeros();
            int[] colOffsets = { 0, nnz };
            int status = SuiteSparsePInvokes.RowAdd(Order, factorizedMatrix, rowIdx,
                nnz, newRow.RawValues, newRow.RawIndices, colOffsets, common);
            if (status != 1)
            {
                throw new SuiteSparseException("Rows addition did not succeed. This could be caused by insufficent memory");
            }
        }

        /// <summary>
        /// Solves the linear system L^T * x = b (or  D * L^T * x = b), where L is the lower triangular factor (and D the 
        /// diagonal factor) of the Cholesky factorization: A = L * L^T (or A = L * D * L^T).
        /// </summary>
        /// <param name="rhsVector">The right hand side vector b of the linear system. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == this.<see cref="Order"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rhsVector"/> violates the described 
        ///     constraints.</exception>
        /// <exception cref="AccessViolationException">Thrown if the unmanaged memory that holds the factorization data has been 
        ///     released.</exception>
        /// <exception cref="SuiteSparseException">Thrown if the call to SuiteSparse library fails.</exception>
        public Vector BackSubstitution(Vector rhsVector)
        {
            var solution = new double[rhsVector.Length];
            SolveInternal(SystemType.BackSubstitution, rhsVector, solution);
            return Vector.CreateFromArray(solution);
        }

    /// <summary>
    /// Solves a series of linear systems L^T * x = b (or  D * L^T * x = b), where L is the lower triangular factor (and D  
    /// the diagonal factor) of the Cholesky factorization: A = L * L^T (or A = L * D * L^T).
    /// </summary>
    /// <param name="rhsVectors">A matrix whose columns are the right hand side vectors b of the linear systems. Constraints:
    ///     <paramref name="rhsVectors"/>.<see cref="IIndexable2D.NumRows"/> == this.<see cref="Order"/>.</param>
    /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rhsVectors"/> violates the described 
    ///     constraints.</exception>
    /// <exception cref="AccessViolationException">Thrown if the unmanaged memory that holds the factorization data has been 
    ///     released.</exception>
    /// <exception cref="SuiteSparseException">Thrown if the call to SuiteSparse library fails.</exception>
    public Matrix BackSubstitutions(Matrix rhsVectors) => SolveInternal(SystemType.BackSubstitution, rhsVectors);

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        /// <exception cref="AccessViolationException">Thrown if the unmanaged memory that holds the factorization data has been 
        ///     released.</exception>
        public double CalcDeterminant()
        {
            if (factorizedMatrix == IntPtr.Zero)
            {
                throw new AccessViolationException("The factorized matrix has been freed from unmanaged memory");
            }
            throw new NotImplementedException();
        }

        /// <summary>
        /// Sets the <paramref name="rowIdx"/>-th row/column of the factorized matrix to the one it would have if the
        /// <paramref name="rowIdx"/>-th row/column of the identity matrix was set as the <paramref name="rowIdx"/>-th 
        /// row/column of the original matrix and then the factorization was performed.
        /// </summary>
        /// <param name="rowIdx">The index of the row/column to modify. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; this.<see cref="Order"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowIdx"/> violates the described constraints.
        ///     </exception>
        /// <exception cref="SuiteSparseException">Thrown if the call to SuiteSparse library fails, usually due to insufficient 
        ///     memory.</exception>
        public void DeleteRow(int rowIdx)
        {
            //TODO: use Preconditions for these tests and implement IIndexable2D.
            if ((rowIdx < 0) || (rowIdx >= Order))
            {
                throw new IndexOutOfRangeException($"Cannot access row {rowIdx} in a"
                    + $" {Order}-by-{Order} matrix");
            }

            int status = SuiteSparsePInvokes.RowDelete(factorizedMatrix, rowIdx, common);
            if (status != 1)
            {
                throw new SuiteSparseException("Rows deletion did not succeed.");
            }
        }

        /// <summary>
        /// See <see cref="IDisposable.Dispose"/>.
        /// </summary>
        public void Dispose()
        {
            ReleaseResources();
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// Solves the linear system L * x = b (or L * D * x = b), where L is the lower triangular factor (and D the 
        /// diagonal factor) of the Cholesky factorization: A = L * L^T (or A = L * D * L^T).
        /// </summary>
        /// <param name="rhsVector">The right hand side vector b of the linear system. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == this.<see cref="Order"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rhsVector"/> violates the described 
        ///     constraints.</exception>
        /// <exception cref="AccessViolationException">Thrown if the unmanaged memory that holds the factorization data has been 
        ///     released.</exception>
        /// <exception cref="SuiteSparseException">Thrown if the call to SuiteSparse library fails.</exception>
        public Vector ForwardSubstitution(Vector rhsVector)
        {
            var solution = new double[rhsVector.Length];
            SolveInternal(SystemType.ForwardSubstitution, rhsVector, solution);
            return Vector.CreateFromArray(solution);
        }

        /// <summary>
        /// Solves a series of linear systems L * x = b (or L * D * x = b), where L is the lower triangular factor (and D  
        /// the diagonal factor) of the Cholesky factorization: A = L * L^T (or A = L * D * L^T).
        /// </summary>
        /// <param name="rhsVectors">A matrix whose columns are the right hand side vectors b of the linear systems. Constraints:
        ///     <paramref name="rhsVectors"/>.<see cref="IIndexable2D.NumRows"/> == this.<see cref="Order"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rhsVectors"/> violates the described 
        ///     constraints.</exception>
        /// <exception cref="AccessViolationException">Thrown if the unmanaged memory that holds the factorization data has been 
        ///     released.</exception>
        /// <exception cref="SuiteSparseException">Thrown if the call to SuiteSparse library fails.</exception>
        public Matrix ForwardSubstitutions(Matrix rhsVectors) => SolveInternal(SystemType.ForwardSubstitution, rhsVectors);

        /// <summary>
        /// See <see cref="ITriangulation.SolveLinearSystem(Vector, Vector)"/>.
        /// </summary>
        /// <remarks>
        /// Solves the linear system L * L^T * x = b (or L * D * L^T * x = b), where L is the lower triangular factor (and D the 
        /// diagonal factor) of the Cholesky factorization: A = L * L^T (or A = L * D * L^T).
        /// </remarks>
        /// 
        /// <exception cref="AccessViolationException">Thrown if the unmanaged memory that holds the factorization data has been 
        ///     released.</exception>
        /// <exception cref="SuiteSparseException">Thrown if the call to SuiteSparse library fails.</exception>
        public void SolveLinearSystem(Vector rhsVector, Vector solution)
        {
            Preconditions.CheckMultiplicationDimensions(Order, solution.Length);
            SolveInternal(SystemType.Regular, rhsVector, solution.RawData);
        }

        /// <summary>
        /// Solves a series of linear systems L * L^T * x = b (or L * D * L^T * x = b), where L is the lower triangular factor   
        /// (and D the diagonal factor) of the Cholesky factorization: A = L * L^T (or A = L * D * L^T).
        /// </summary>
        /// <param name="rhsVectors">A matrix whose columns are the right hand side vectors b of the linear systems. Constraints:
        ///     <paramref name="rhsVectors"/>.<see cref="IIndexable2D.NumRows"/> == this.<see cref="Order"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rhsVectors"/> violates the described 
        ///     constraints.</exception>
        /// <exception cref="AccessViolationException">Thrown if the unmanaged memory that holds the factorization data has been 
        ///     released.</exception>
        /// <exception cref="SuiteSparseException">Thrown if the call to SuiteSparse library fails.</exception>
        public Matrix SolveLinearSystems(Matrix rhsVectors) => SolveInternal(SystemType.Regular, rhsVectors);

        /// <summary>
        /// Perhaps I should use SafeHandle (thread safety, etc). 
        /// Also perhaps there should be dedicated objects for closing each handle.
        /// </summary>
        private void ReleaseResources() 
        {
            if (common != IntPtr.Zero)
            {
                // Supposedly throwing in destructors and Dispose() is poor practice.
                if (factorizedMatrix == IntPtr.Zero) 
                {
                    throw new AccessViolationException("The matrix in unmanaged memory has already been cleared or lost");
                }
                SuiteSparsePInvokes.DestroyFactor(ref factorizedMatrix, common);
                factorizedMatrix = IntPtr.Zero;
                SuiteSparsePInvokes.DestroyCommon(ref common);
                common = IntPtr.Zero;
            }
        }

        private void SolveInternal(SystemType system, Vector rhs, double[] solution)
        {
            if (factorizedMatrix == IntPtr.Zero)
            {
                throw new AccessViolationException("The factorized matrix has been freed from unmanaged memory");
            }
            Preconditions.CheckSystemSolutionDimensions(Order, rhs.Length);

            int status = SuiteSparsePInvokes.Solve((int)system, Order, 1, factorizedMatrix, rhs.RawData, solution, common);
            if (status != 1) throw new SuiteSparseException("System solution failed.");
        }

        private Matrix SolveInternal(SystemType system, Matrix rhs)
        {
            if (factorizedMatrix == IntPtr.Zero)
            {
                throw new AccessViolationException("The factorized matrix has been freed from unmanaged memory");
            }
            Preconditions.CheckSystemSolutionDimensions(Order, rhs.NumRows);
            double[] solution = new double[rhs.NumRows * rhs.NumColumns];
            int status = SuiteSparsePInvokes.Solve((int)system, Order, rhs.NumColumns, factorizedMatrix, rhs.RawData,
                solution, common);
            if (status != 1) throw new SuiteSparseException("System solution failed.");
            return Matrix.CreateFromArray(solution, rhs.NumRows, rhs.NumColumns, false);
        }
    }
}
