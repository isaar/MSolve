using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: also return the nonzeros after cholesky, flop count and other statistics
namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    public class OrderingCAMD
    {
        private readonly bool aggressiveAbsorption;
        private readonly int denseThreshold;

        public OrderingCAMD(int denseThreshold = -1, bool aggressiveAbsorption = true)
        {
            this.denseThreshold = denseThreshold;
            this.aggressiveAbsorption = aggressiveAbsorption;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="order"></param>
        /// <param name="nonZerosUpper"></param>
        /// <param name="cscRowIndices"></param>
        /// <param name="cscColOffsets"></param>
        /// <param name="constraints">Array of length = order with ordering constraints. Its values must be 
        ///     0 &lt;= <paramref name="constraints"/>[i] &lt; order. If <paramref name="constraints"/> = NULL, no constraints 
        ///     will be enforced.
        ///		Example: <paramref name="constraints"/> = { 2, 0, 0, 0, 1 }. This means that indices 1, 2, 3 that have 
        ///		<paramref name="constraints"/>[i] = 0, will be ordered before index 4 with <paramref name="constraints"/>[4] = 1,
        ///		which will be ordered before index 0 with <paramref name="constraints"/>[0] = 2. Indeed for a certain pattern, 
        ///		permutation = { 3, 2, 1, 4, 0 } (remember permutation is a new-to-old 
        ///		mapping).</param>
        /// <returns></returns>
        public (int[] permutation, ReorderingStatistics stats) FindPermutation(int order, int nonZerosUpper, int[] cscRowIndices,
            int[] cscColOffsets, int[] constraints)
        {
            var permutation = new int[order];
            int status = SuiteSparseUtilities.ReorderCAMD(order, cscRowIndices, cscColOffsets, constraints, 
                denseThreshold, aggressiveAbsorption ? 1 : 0, 
                permutation, out int nnzFactor, out int numMovedDense);

            // Error checking
            //if (status == 1) throw new SuiteSparseException("The matrix had unsorted columns, but was otherwise ok.");
            if (status == 2) throw new SuiteSparseException(
                "The arrays that describe the non zero pattern of the matrix were invalid or the permutation buffer was NULL.");
            else if (status == 3) throw new SuiteSparseException("Not enough memory could be allocated");
            return (permutation, new ReorderingStatistics(nnzFactor, numMovedDense));
        }
    }
}
