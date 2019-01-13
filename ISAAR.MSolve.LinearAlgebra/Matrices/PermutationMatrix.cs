using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Efficient implementation for permutation matrices, namely matrices that change the order of vector entries, matrix rows
    /// or matrix columns.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PermutationMatrix
    {
        private readonly int order;
        private readonly int[] data;

        private PermutationMatrix(int order, int[] data)
        {
            this.order = order;
            this.data = data;
        }

        /// <summary>
        /// The number of columns of the permutation matrix. 
        /// </summary>
        public int NumColumns { get { return order; } }

        /// <summary>
        /// The number of rows of the permutation matrix. 
        /// </summary>
        public int NumRows { get { return order; } }

        /// <summary>
        /// Initializes a new instance of <see cref="PermutationMatrix"/> that is equal to the identity matrix, namely a square 
        /// matrix with non-diagonal entries being equal to 0 and diagonal entries being equal to 1. This permutation does not
        /// change modify the vector/matrix it is applied to. However, it can be used as a starting point to define other 
        /// permutations.
        /// </summary>
        /// <param name="order">The number of rows/columns of the identity matrix.</param>
        public static PermutationMatrix CreateIdentity(int order)
        {
            var data = new int[order];
            for (int i = 0; i < order; ++i) data[i] = i;
            return new PermutationMatrix(order, data);
        }

        /// <summary>
        /// Modifies this <see cref="PermutationMatrix"/> instance, such that it defines a row-exchange operation: 
        /// <paramref name="rowIdx1"/> becomes <paramref name="rowIdx2"/> and vice versa.
        /// </summary>
        /// <param name="rowIdx1">The index of the first row to exchange.</param>
        /// <param name="rowIdx2">The index of the second row to exchange.</param>
        public void ExchangeRows(int rowIdx1, int rowIdx2)
        {
            int swap = data[rowIdx1];
            data[rowIdx1] = data[rowIdx2];
            data[rowIdx2] = swap; 
        }

        /// <summary>
        /// Multiplies this permutation matrix or its transpose with a vector: result = oper(this) * <paramref name="vector"/>.
        /// This is equivalent to applying the permutation defined by this <see cref="PermutationMatrix"/> to 
        /// <paramref name="vector"/>. If (<paramref name="transposeThis"/> == true) this permutation is new-to-old:
        /// result[i] = <paramref name="vector"/>[permutation[i]]. Otherwise it is old-to-new:
        /// result[permutation[i]] = <paramref name="vector"/>[i]
        /// </summary>
        /// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to the <see cref="NumColumns"/>
        ///     of oper(this).</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
        ///     <paramref name="vector"/> is different than the <see cref="NumColumns"/> of oper(this).</exception>
        public Vector MultiplyRight(Vector vector, bool transposeThis)
        {
            Preconditions.CheckMultiplicationDimensions(order, vector.Length);
            var result = new double[order];
            if (transposeThis)
            {
                for (int i = 0; i < order; ++i) result[i] = vector[data[i]]; // Verify this
            }
            else
            {
                for (int i = 0; i < order; ++i) result[data[i]] = vector[i];
            }
            return Vector.CreateFromArray(result, false);
        }

    }
}
