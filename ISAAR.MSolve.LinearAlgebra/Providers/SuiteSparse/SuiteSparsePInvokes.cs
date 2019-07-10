using System;
using System.Runtime.InteropServices;

//TODO: Use enums as return values or at least named constants.
namespace ISAAR.MSolve.LinearAlgebra.Providers.PInvoke
{
    /// <summary>
    /// Platform Invoke methods for calling the SuiteSparse library via my custom C interface.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SuiteSparsePInvokes
    {
        /// <summary>
        /// Allocates in heap and returns a handle with matrix settings that must be passed to all CHOLMOD functions. Returns
        /// <see cref="IntPtr.Zero"/> if any failure occurs.
        /// </summary>
        /// <param name="factorizationType">0 for for simplicial L^T*L or L^T*D*L factorization, 
        ///     1 for supernodal L^T*L factorization, 
        ///     2 for automatic decision between supernodal/simplicial factorization, 
        ///     3 for automatic decision and after factorization convert to simplicial.
        ///     Supernodal is usually faster and results in faster solutions, but to modify the factorized matrix it must be
        ///     converted to simplicial, though this can be	done automatically.</param>
        /// <param name="orderingType">0 for no reordering, 1 for automatic reordering (let suitesparse try some alternatives and
        ///     keep the best), 2 for AMD.</param>
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
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_factorize_cscupper")]
        internal static extern int FactorizeCSCUpper(int order, int nnz, double[] values, int[] rowIndices, int[] colOffsets, 
            out IntPtr factorizedMatrix, IntPtr common);

        /// <summary>
        /// Returns the number of non zero entries in the factorized matrix (only the explitily stored factor). If anything goes 
        /// wrong -1 is returned.
        /// </summary>
        /// <param name="factorization">Pointer to the factorized matrix, which is stored in unmanaged memory.</param>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_get_factor_nonzeros")]
        internal static extern int GetFactorNonZeros(IntPtr factorization);

        /// <summary>
        /// Caclulates a fill reducing ordering using the Approximate Minimum Degree algorithm for a symmetric sparse matrix.  
        /// Returns 1 if the reordering is successful, 0 if it failed (e.g. due to exceeding the available memory).
        /// </summary> 
        /// <param name="order">Number of rows = number of columns.</param>
        /// <param name="nnz">Number of non zero entries in the upper triangle.</param>
        /// <param name="rowIndices">Array containing the row indices of the non zero entries of the upper triangle. 
        ///     Length = <paramref name="nnz"/>. They must be sorted.</param>
        /// <param name="colOffsets">Array containing the indices into <paramref name="rowIndices"/> of the first entry of each 
        ///     column. Length = <paramref name="order"/> + 1. They must be sorted. 
        ///     The first entry is <paramref name="colOffsets"/>[0] = 0. 
        ///     The last entry is <paramref name="colOffsets"/>[<paramref name="order"/>] = nnz.</param>
        /// <param name="outPermutation">Out parameter: buffer of length = <paramref name="order"/>. Will be filled with a  
        ///     fill-reducing permutation vector, such that: original index = i, reordered index = 
        ///     <paramref name="outPermutation"/>[i].</param>
        /// <param name="outFactorNNZ">Out parameter: the number of non zero entries in a subsequent L*L^T factorization. Will 
        ///     be -1 if the ordering fails.</param>
        /// <param name="common">The matrix settings.</param>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_reorder_amd_upper")]
        internal static extern int ReorderAMDUpper(int order, int nnz, int[] rowIndices, int[] colOffsets, int[] outPermutation,
            [Out] out int outFactorNNZ, IntPtr common);

        /// <summary>
        /// Calculates a fill reducing ordering using the Constrained Approximate Minimum Degree algorithm for a A + A^T,  
        /// where A is a square sparse matrix. The pattern of A + A^T is formed first. The constrains enforce groups of indices 
        /// to be ordered consecutively, before other groups.
        /// Returns: 
        ///	    0 if the input was ok and the ordering is successful, 
        ///	    1 if the matrix had unsorted columns or duplicate entries, but was otherwise valid, 
        ///	    2 if input arguments <paramref name="order"/>, <paramref name="colOffsets"/>, <paramref name="rowIndices"/> are 
        ///	        invalid, or if <paramref name="outPermutation"/> is NULL,
        ///	    3 if not enough memory can be allocated.
        /// </summary> 
        /// <param name="order">Number of rows = number of columns.</param>
        /// <param name="rowIndices">Array containing the row indices of the non zero entries of the upper triangle. 
        ///     Length = <paramref name="nnz"/>. They must be sorted.</param>
        /// <param name="colOffsets">Array containing the indices into <paramref name="rowIndices"/> of the first entry of each 
        ///     column. Length = <paramref name="order"/> + 1. They must be sorted. 
        ///     The first entry is <paramref name="colOffsets"/>[0] = 0. 
        ///     The last entry is <paramref name="colOffsets"/>[<paramref name="order"/>] = nnz.</param>
        /// <param name="constraints">Array of length = order with ordering constraints. Its values must be 
        ///     0 &lt;= <paramref name="constraints"/>[i] &lt; order. If <paramref name="constraints"/> = NULL, no constraints 
        ///     will be enforced.
        ///		Example: <paramref name="constraints"/> = { 2, 0, 0, 0, 1 }. This means that indices 1, 2, 3 that have 
        ///		<paramref name="constraints"/>[i] = 0, will be ordered before index 4 with <paramref name="constraints"/>[4] = 1,
        ///		which will be ordered before index 0 with <paramref name="constraints"/>[0] = 2. Indeed for a certain pattern, 
        ///		<paramref name="outPermutation"/> = { 3, 2, 1, 4, 0 } (remember <paramref name="outPermutation"/> is a new-to-old 
        ///		mapping).</param>
        ///	<param name="denseThreshold">A dense row/column in A + A^T can cause CAMD to spend significant time in ordering    
        /// 	the matrix. If <paramref name="denseThreshold"/> &gt;= 0, rows/columns with more than 
        /// 	<paramref name="denseThreshold"/> * sqrt(order(A)) entries are ignored during the ordering, and placed last in  
        /// 	the output order. The default value of <paramref name="denseThreshold"/> is 10. If negative, no rows/columns are 
        /// 	treated as dense. Rows/columns with 16 or fewer off-diagonal entries are never considered dense. WARNING: 
        /// 	allowing dense rows/columns may violate the constraints.</param>
        /// <param name="aggressiveAbsorption">If non zero, aggressive absorption will be performed, which means that a  
        /// 	prior element is absorbed into the current element if it is a subset of the current element, even if it is not  
        /// 	adjacent to the current pivot element. This nearly always leads to a better ordering (because the approximate  
        /// 	degrees are more accurate) and a lower execution time. There are cases where it can lead to a slightly worse 
        /// 	ordering, however. The default value is nonzero. To turn it off, set <paramref name="aggressiveAbsorption"/> 
        /// 	to 0.</param>
        /// <param name="outPermutation">Out parameter: buffer of length = <paramref name="order"/>. Will be filled with a  
        ///     fill-reducing permutation vector, such that: original index = i, reordered index = 
        ///     <paramref name="outPermutation"/>[i].</param>
        /// <param name="outFactorNNZ">Out parameter: upper bound on the number of non zero entries in L of a subsequent L*L^T   
        ///     factorization. Will be -1 if the ordering fails.</param>
        /// <param name="outMovedDense">Out parameter: the number of dense rows/columns of A + A^T that were removed from A  
        ///     prior to ordering. These are placed last in the output order of <paramref name="outPermutation"/>. Will be -1 
        ///     if the ordering fails. WARNING: if <paramref name="outMovedDense"/> &gt; 0, it indicates that the constraints  
        ///     are violated.</param>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_reorder_camd")]
        internal static extern int ReorderCAMD(int order, int[] rowIndices, int[] colOffsets, int[] constraints,
            int denseThreshold, int aggressiveAbsorption, int[] outPermutation, [Out] out int outFactorNNZ, 
            [Out] out int outMovedDense);

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
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_row_delete")]
        internal static extern int RowDelete(IntPtr factorizedMatrix, int rowIdx, IntPtr common);

        /// <summary>
        /// Solves a linear system or applies back substitution or forward substituiton to 1 ore more right hand sides.
        /// Returns 1 if the method succeeds, 0 otherwise.
        /// </summary>
        /// <param name="system">0 for system solution (A*x=b), 4 for forward substitution (L*x=b), 
        ///     5 for back substitution (L^T*x=b)</param>
        /// <param name="numRows">Number of matrix rows = number of matrix columns = number of rhs matrix rows.</param>
        /// <param name="numRhs">Number of rhs vectors = number of columns in rhs matrix.</param>
        /// <param name="factorizedMatrix">The data of the cholesky factorization of the matrix.</param>
        /// <param name="rhs">The right hand side matrix. Column major array with dimensions = 
        ///     <paramref name="numRows"/> -by- <paramref name="numRhs"/>.</param>
        /// <param name="outSolution">Buffer for the left hand side vector (unknown). Column major array with dimensions = 
        /// 	<paramref name="numRows"/> -by- <paramref name="numRhs"/>.</param>
        /// <param name="common">The matrix settings.</param>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_solve")]
        internal static extern int Solve(int system, int numRows, int numRhs, IntPtr factorizedMatrix, 
            double[] rhs, double[] outSolution, IntPtr common);
    }
}