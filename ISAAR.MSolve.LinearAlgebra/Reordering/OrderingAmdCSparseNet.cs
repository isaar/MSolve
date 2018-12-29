using System;
using CSparse;
using CSparse.Double;
using CSparse.Ordering;

//TODO: Creating a dummy variables array with as many entries as the non zero entries is expensive and not needed. I would be 
//      better off copying the AMD source code and avoiding that step.
//TODO: Find out what is going wrong when AMD returns a permutation with more entries than the matrix order.
namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    /// <summary>
    /// Calculates a fill-reducing permutation for the rows/columns of a symmetric sparse matrix, using the Approximate Minimum 
    /// Degree (AMD) ordering algorithm. The AMD implementation used is provided by the CSparse.NET library. For more 
    /// information, see the AMD user guide, which is distributed as part of the SuiteSparse library.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class OrderingAmdCSparseNet : IReorderingAlgorithm
    {
        /// <summary>
        /// See <see cref="IReorderingAlgorithm.FindPermutation(SparsityPatternSymmetric)"/>
        /// </summary>
        /// <remarks>The returned permutation is new-to-old.</remarks>
        public (int[] permutation, bool oldToNew) FindPermutation(SparsityPatternSymmetric pattern)
        {
            int order = pattern.Order;
            (int[] cscRowIndices, int[] cscColOffsets) = pattern.BuildSymmetricCSCArrays(true); //TODO: perhaps sorting is not needed here.
            var dummyCscValues = new double[cscRowIndices.Length]; //TODO: too expensive 
            var matrixCSparse = new SparseMatrix(order, order, dummyCscValues, cscRowIndices, cscColOffsets);
            int[] permutation = AMD.Generate<double>(matrixCSparse, ColumnOrdering.MinimumDegreeAtPlusA);

            // It is possible that CSparse.NET AMD algorithm returns more entries than the matrix order (so far I have found 1 
            // extra). In that case, make sure the first ones are valid and return only them.
            if (permutation.Length > order)
            {
                for (int i = order; i < permutation.Length; ++i)
                {
                    if (permutation[i] < pattern.Order) throw new Exception(
                        "Something went wrong during AMD. The permutation vector has more entries than the matrix order.");
                }
                var permutationCorrected = new int[order];
                Array.Copy(permutation, permutationCorrected, order);
                return (permutationCorrected, false);
            }
            else return (permutation, false);
        }
    }
}
