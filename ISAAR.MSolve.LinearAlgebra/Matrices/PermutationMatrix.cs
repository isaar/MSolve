using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public class PermutationMatrix
    {
        private readonly int order;
        private readonly int[] data;

        private PermutationMatrix(int order, int[] data)
        {
            this.order = order;
            this.data = data;
        }

        public int NumColumns { get { return order; } }
        public int NumRows { get { return order; } }

        public static PermutationMatrix CreateIdentity(int order)
        {
            var data = new int[order];
            for (int i = 0; i < order; ++i) data[i] = i;
            return new PermutationMatrix(order, data);
        }

        public void ExchangeRows(int idx1, int idx2)
        {
            int swap = data[idx1];
            data[idx1] = data[idx2];
            data[idx2] = swap; 
        }

        public Vector MultiplyRight(Vector vector, bool transpose)
        {
            Preconditions.CheckMultiplicationDimensions(order, vector.Length);
            var result = new double[order];
            if (transpose)
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
