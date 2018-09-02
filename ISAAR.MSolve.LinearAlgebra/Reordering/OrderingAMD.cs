using System;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;

//TODO: also return the nonzeros after cholesky, flop count and other statistics
//TODO: the number of "dense" rows moved to the end is not reported by SuiteSparse and a dummy value (-1) is return. Fix this.
namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    /// <summary>
    /// Implements the Approximate Minimum Degree (AMD) ordering algorithm, by calling the appropriate SuiteSparse library 
    /// functions. For more, see the AMD user guide, which is distributed as part of the SuiteSparse library.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class OrderingAmd
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="OrderingAmd"/> class.
        /// </summary>
        public OrderingAmd()
        {
            //TODO: add options for dense rows and aggressive absorption.
        }

        /// <summary>
        /// Find a fill reducing permutation for the sparsity pattern of a symmetric matrix defined by the parameters.
        /// The returned permutation is new-to-old, namely reordered[i] = original[permutation[i]].
        /// </summary>
        /// <param name="order">The number of rows/columns of the symmetric matrix.</param>
        /// <param name="nonZerosUpper">The number of (structural) non-zero entries in the upper triangle of the symmetric 
        ///     matrix.</param>
        /// <param name="cscRowIndices">Row indices of the upper triangle entries of the symmetric matrix, in Compressed Sparse 
        ///     Columns format. All row indices of the same column must be sorted.</param>
        /// <param name="cscColOffsets">Column offsets of the upper triangle entries of the symmetric matrix, in Compressed  
        ///     Sparse Columns format. All column offsets must be sorted.</param>
        /// <returns>permutation: An array containing the new-to-old fill reducing permutation. 
        ///     stats: Measuments taken by SuiteSparse during the execution of AMD.</returns>
        /// <exception cref="SuiteSparseException">Thrown if SuiteSparse dlls cannot be loaded or if AMD fails.</exception>
        public (int[] permutation, ReorderingStatistics stats) FindPermutation(int order, int nonZerosUpper, int[] cscRowIndices, 
            int[] cscColOffsets)
        {
            var permutation = new int[order];
            IntPtr common = SuiteSparseUtilities.CreateCommon(0, 0);
            if (common == IntPtr.Zero) throw new SuiteSparseException("Failed to initialize SuiteSparse.");
            int status = SuiteSparseUtilities.ReorderAMDUpper(order, nonZerosUpper, cscRowIndices, cscColOffsets, permutation, 
                out int nnzFactor, common);
            if (status == 0) throw new SuiteSparseException("AMD failed. This could be caused by the matrix being so large it"
                + " cannot be processed with the available memory.");
            SuiteSparseUtilities.DestroyCommon(ref common);
            return (permutation, new ReorderingStatistics(nnzFactor, -1));
        }

        /// <summary>
        /// Find a fill reducing permutation for the sparsity pattern of a symmetric matrix.
        /// The returned permutation is new-to-old, namely reordered[i] = original[permutation[i]].
        /// </summary>
        /// <param name="pattern">The indices of the non-zero entries of a symmetric matrix.</param>
        /// <returns>permutation: An array containing the new-to-old fill reducing permutation. 
        ///     stats: Measuments taken by SuiteSparse during the execution of AMD.</returns>
        /// <exception cref="SuiteSparseException">Thrown if SuiteSparse dlls cannot be loaded or if AMD fails.</exception>
        public (int[] permutation, ReorderingStatistics stats) FindPermutation(SparsityPatternSymmetric pattern)
        {
            (int[] rowIndices, int[] colOffsets) = pattern.BuildSymmetricCSCArrays(true);
            return FindPermutation(pattern.Order, rowIndices.Length, rowIndices, colOffsets);
        }

        /// <summary>
        /// Find a fill reducing permutation for the sparsity pattern of a symmetric matrix.
        /// The returned permutation is new-to-old, namely reordered[i] = original[permutation[i]].
        /// </summary>
        /// <param name="dok">A builder for symmetric sparse matrices.</param>
        /// <returns>permutation: An array containing the new-to-old fill reducing permutation. 
        ///     stats: Measuments taken by SuiteSparse during the execution of AMD.</returns>
        /// <exception cref="SuiteSparseException">Thrown if SuiteSparse dlls cannot be loaded or if AMD fails.</exception>
        public (int[] permutation, ReorderingStatistics stats) FindPermutation(DokSymmetric dok)
        {
            (double[] values, int[] rowIndices, int[] colOffsets) = dok.BuildSymmetricCscArrays(true);
            return FindPermutation(dok.NumColumns, rowIndices.Length, rowIndices, colOffsets);
        }
    }
}
