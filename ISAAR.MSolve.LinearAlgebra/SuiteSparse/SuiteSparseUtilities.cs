using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.SuiteSparse
{
    internal static class SuiteSparseUtilities
    {
        /// <summary>
        /// Allocates in heap and returns a handle with matrix settings that must be passed to all CHOLMOD functions. Returns
        /// <see cref="IntPtr.Zero"/> if any failure occurs.
        /// </summary>
        /// <param name="factorizationType">0 for for simplicial factorization, 1 for automatic decidion between supernodal / 
        ///     simplicial factorization, 2 for automatic decidion between supernodal/simplicial factorization and after 
        ///     factorization convert to simplicial. Supernodal is usually faster, but to modify the factorized matrix it must be
        ///     converted to simplicial, though this can be	done automatically.</param>
        /// <param name="orderingType">0 for no reordering, 1 for automatic reordering (let suitesparse try some alternatives and
        ///     keep the best).</param>
        /// <returns></returns>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_create_common")]
        internal static extern IntPtr CreateCommon(int factorizationType, int orderingType);

        /// <summary>
        /// Frees the memory for the matrix settings in unmanaged memory.
        /// </summary>
        /// <param name="common">The matrix settings. It will be freed in unmanaged memory.</param>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_destroy_common")]
        internal static extern void DestroyCommon(ref IntPtr common);

        /// <summary>
        /// Frees the memory for the factorized matrix in unmanaged memory
        /// </summary>
        /// <param name="factorizedMatrix">The factorized matrix data. It will be freed in unmanaged memory.</param>
        /// <param name="common">The matrix settings.</param>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_destroy_factor")]
        internal static extern void DestroyFactor(ref IntPtr factorizedMatrix, IntPtr common);

        /// <summary>
        /// Factorize a symmetric matrix using cholesky algorithm. The matrix is in csc form, with only the upper triangle stored.
        /// If cholesky is successful: -1 is returned and factorizedMatrix points to the factorized upper triangle.
        /// If the matrix is not positive definite: the index(0-based) of the column where cholesky failed
        /// and factorizedMatrix = <see cref="IntPtr.Zero"/>.
        /// If the something another failure occurs, such as memory not being sufficient due to excessive fill-in: -2 is returned
        /// and factorizedMatrix = <see cref="IntPtr.Zero"/>.
        /// </summary>
        /// <param name="order">Number of rows = number of columns.</param>
        /// <param name="nnz">Number of non zero entries in the upper triangle.</param>
        /// <param name="values">Array containing the non zero entries of the upper triangle in column major order. 
        ///     Length = <paramref name="nnz"/>.</param>
        /// <param name="rowIndices">Array containing the row indices of the non zero entries of the upper triangle. 
        ///     <paramref name="nnz"/>. They must be sorted.</param>
        /// <param name="colOffsets">Array containing the indices into <paramref name="values"/> (and 
        ///     <paramref name="rowIndices"/>) of the first entry of each column. Length = <paramref name="order"/> + 1. They 
        ///     must be sorted. The first entry is <paramref name="colOffsets"/>[0] = 0. The last entry is 
        ///     <paramref name="colOffsets"/>[<paramref name="order"/>] = <paramref name="nnz"/>.</param>
        /// <param name="factorizedMatrix"> Out parameter - the factorized upper triangle of the symmetric CSC matrix.</param>
        /// <param name="common">The matrix settings.</param>
        /// <returns></returns>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_factorize_cscupper")]
        internal static extern int FactorizeCSCUpper(int order, int nnz, double[] values, int[] rowIndices, int[] colOffsets, 
            out IntPtr factorizedMatrix, IntPtr common);

        /// <summary>
        /// Returns the number of non zero entries in the factorized matrix. If anything goes wrong -1 is returned.
        /// </summary>
        /// <param name="factorization">Pointer to the factorized matrix, which is stored in unmanaged memory.</param>
        /// <returns></returns>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_get_factor_nonzeros")]
        internal static extern int GetFactorNonZeros(IntPtr factorization);

        /// <summary>
        /// 
        /// </summary> 
        /// <param name="order">Number of rows = number of columns.</param>
        /// <param name="nnz">Number of non zero entries in the upper triangle.</param>
        /// <param name="rowIndices">Array containing the row indices of the non zero entries of the upper triangle. 
        ///     Length = <paramref name="nnz"/>. They must be sorted.</param>
        /// <param name="colOffsets">Array containing the indices into <paramref name="rowIndices"/> of the first entry of each 
        ///     column. Length = <paramref name="order"/> + 1. They must be sorted. 
        ///     The first entry is <paramref name="colOffsets"/>[0] = 0. 
        ///     The last entry is <paramref name="colOffsets"/>[<paramref name="order"/>] = nnz.</param>
        /// <param name="outPermutation">Buffer of length = <paramref name="order"/>. Will be filled with a fill-reducing 
        ///     permutation vector, such that: original index = i, reordered index = <paramref name="outPermutation"/>[i].</param>
        /// <param name="common">The matrix settings.</param>
        /// <returns></returns>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_reorder_amd_upper")]
        internal static extern int ReorderAMDUpper(int order, int nnz, int[] rowIndices, int[] colOffsets, int[] outPermutation,
            IntPtr common);

        /// <summary>
        /// Adds a row and column to an LDL' factorization. Before updating the kth row and column of L must be equal to the kth  
        /// row and column of the identity matrix. The row/column to add must be a sparse CSC matrix with dimensions n-by-1,  
        /// where n is the order of the matrix. Returns 1 if the method succeeds, 0 otherwise.
        /// </summary>
        /// <param name="order">Number of rows = number of columns.</param>
        /// <param name="factorizedMatrix">The data of the cholesky factorization of the matrix. It will be modified.</param>
        /// <param name="rowIdx">Index of row/column to add.</param>
        /// <param name="vectorNnz">Number of non zero entries of the row/column to add.</param>
        /// <param name="vectorValues">The CSC format values of the row/column to add.</param>
        /// <param name="vectorRowIndices">The CSC format row indices of the row/column to add.</param>
        /// <param name="vectorColOffsets">The CSC format column offsets of the row/column to add.</param>
        /// <param name="common">The matrix settings.</param>
        /// <returns></returns>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_row_add")]
        internal static extern int RowAdd(int order, IntPtr factorizedMatrix, int rowIdx, int vectorNnz, 
            double[] vectorValues, int[] vectorRowIndices, int[] vectorColOffsets, IntPtr common);

        /// <summary>
        /// Deletes a row and column from an cholesky factorization. After updating the kth row and column of L will be equal to 
        /// the kth row and column of the identity matrix.Returns 1 if the method succeeds, 0 otherwise.
        /// </summary>
        /// <param name="factorizedMatrix">The LDL' factorization of the matrix. It will be modified.</param>
        /// <param name="rowIdx">Index of row/column to delete.</param>
        /// <param name="common">The matrix settings.</param>
        /// <returns></returns>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_row_delete")]
        internal static extern int RowDelete(IntPtr factorizedMatrix, int rowIdx, IntPtr common);

        /// <summary>
        /// Solves a linear system with a single right hand side vector.
        /// </summary>
        /// <param name="order">Number of matrix rows = number of matrix columns = length of right hand side vector.</param>
        /// <param name="factorizedMatrix">The data of the cholesky factorization of the matrix.</param>
        /// <param name="rhs">The right hand side vector. Its length must be equal to the order of the matrix: 
        ///     factorized_matrix->n.</param>
        /// <param name="outSolution">Buffer for the left hand side vector (unknown). Its length must be equal to the order of 
        ///     the matrix: factorized_matrix->n.</param>
        /// <param name="common">The matrix settings.</param>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_solve")]
        internal static extern void Solve(int order, IntPtr factorizedMatrix, double[] rhs, double[] outSolution, IntPtr common);
    }
}