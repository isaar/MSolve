using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.SuiteSparse;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    public class CholeskySuiteSparse : IFactorization, IDisposable
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
        /// The number of rows or columns of the matrix. 
        /// </summary>
        public int Order { get; }

        public static CholeskySuiteSparse Factorize(int order, int nonZerosUpper, double[] values, int[] rowIndices,
            int[] colOffsets)
        {
            IntPtr common = SuiteSparseUtilities.CreateCommon();
            int status = SuiteSparseUtilities.FactorizeCSCUpper(order, nonZerosUpper, values, rowIndices, colOffsets, 
                out IntPtr factorizedMatrix, common);
            if (status == -2)
            {
                SuiteSparseUtilities.DestroyCommon(ref common);
                throw new SuiteSparseException("Factorization did not succeed. This could be caused by insufficent memory,"
                    + " due to excessive fill-in.");
            }
            else if (status >= 0)
            {
                SuiteSparseUtilities.DestroyCommon(ref common);
                throw new SuiteSparseException("The matrix not being positive definite."
                    + $" Cholesky failed at column {status} (0-based indexing).");
            }
            else return new CholeskySuiteSparse(order, common, factorizedMatrix);
        }

        /// <summary>
        /// Update row (same as column) <paramref name="rowIdx"/> of the factorized matrix to the one it would have if 
        /// <paramref name="newRow"/> was set as the <paramref name="rowIdx"/>-th row/col of the original matrix and then the 
        /// factorization was computed. The existing <paramref name="rowIdx"/>-th row/column of the original matrix must be equal 
        /// to the <paramref name="rowIdx"/>-th row/col of the identity matrix. 
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="newRow"></param>
        public void AddRow(int rowIdx, SparseVector newRow) //TODO: The row should be input as a sparse CSC matrix with dimensions order-by-1
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
            int status = SuiteSparseUtilities.RowAdd(Order, factorizedMatrix, rowIdx, 
                nnz, newRow.InternalValues, newRow.InternalRowIndices, colOffsets, common);
            if (status != 1)
            {
                throw new SuiteSparseException("Rows addition did not succeed. This could be caused by insufficent memory");
            }
        }

        public double CalcDeterminant()
        {
            if (factorizedMatrix == IntPtr.Zero)
            {
                throw new AccessViolationException("The factorized matrix has been freed from unmanaged memory");
            }
            throw new NotImplementedException();
        }

        /// <summary>
        /// Update row (same as column) <paramref name="rowIdx"/> of the factorized matrix to the one it would have if the
        /// <paramref name="rowIdx"/>-th row/col of the identity matrix was set as the <paramref name="rowIdx"/>-th row/col of  
        /// the original matrix and then the factorization was computed.
        /// </summary>
        /// <param name="rowIdx"></param>
        public void DeleteRow(int rowIdx)
        {
            int status = SuiteSparseUtilities.RowDelete(factorizedMatrix, rowIdx, common);
            if (status != 1)
            {
                throw new SuiteSparseException("Rows deletion did not succeed.");
            }
        }

        public void Dispose()
        {
            ReleaseResources();
            GC.SuppressFinalize(this);
        }

        public VectorMKL SolveLinearSystem(VectorMKL rhs)
        {
            if (factorizedMatrix == IntPtr.Zero)
            {
                throw new AccessViolationException("The factorized matrix has been freed from unmanaged memory");
            }
            double[] solution = new double[rhs.Length];
            SuiteSparseUtilities.Solve(Order, factorizedMatrix, rhs.InternalData, solution, common);
            return VectorMKL.CreateFromArray(solution, false);
        }

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
                SuiteSparseUtilities.DestroyFactor(ref factorizedMatrix, common);
                factorizedMatrix = IntPtr.Zero;
                SuiteSparseUtilities.DestroyCommon(ref common);
                common = IntPtr.Zero;
            }
        }
    }
}
