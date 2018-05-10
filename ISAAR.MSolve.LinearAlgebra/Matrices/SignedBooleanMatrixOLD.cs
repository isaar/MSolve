using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: this is not correct. 1) if the boolean matrix is for the whole domain, then each row will have exactly 2 nonzero entries
//      a 1 and a -1. 2) Whether the boolean matrix is for the whole domain or a subdomain, columns that correspond to dofs 
//      with multiplicity = 2 will have exactly 1 nonzero entry, while columns that correspond to dofs with multiplicity > 2
//      will have exactly 2 nonzero entries.
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Each row corresponds to a displacement continuity equation. The columns are the dofs of the subdomain. Each row can have
    /// at most one non zero entry (1 or -1) at some column that corresponds to a boundary dof. Also each column can have at most
    /// one non zero entry.
    /// </summary>
    public class SignedBooleanMatrixOLD: IIndexable2D
    {
        private readonly Dictionary<int, (int col, int sign)> data;

        public SignedBooleanMatrixOLD(int numRows, int numColumns)
        {
            this.NumRows = numRows;
            this.NumColumns = numColumns;
            this.data = new Dictionary<int, (int col, int sign)>();
        }

        public int NumRows { get; }
        public int NumColumns { get; }

        public double this[int rowIdx, int colIdx]
        {
            get
            {
                bool rowExists = data.TryGetValue(rowIdx, out (int, int) colSignPair);
                if (rowExists && colIdx == colSignPair.Item1) return colSignPair.Item2;
                return 0.0;
            }
        }

        /// <summary>
        /// Sets B[<paramref name="rowIdx"/>, <paramref name="colIdx"/>] = 1 or -1.
        /// </summary>
        /// <param name="rowIdx">The row index. This row must be empty before a non zero entry is added.</param>
        /// <param name="colIdx">The column index for the row <paramref name="rowIdx"/>.</param>
        /// <param name="sign">True for +1, false for -1</param>
        public void AddEntry(int rowIdx, int colIdx, bool sign)
        {
            if (sign) data.Add(rowIdx, (colIdx, 1));
            else data.Add(rowIdx, (colIdx, -1));
        }

        public double[,] CopyToArray2D()
        {
            return DenseStrategies.CopyToArray2D(this);
        }

        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        public Vector MultiplyRight(Vector vector, bool transposeThis)
        {
            if (transposeThis) return MultiplyTransposed(vector);
            else return MultiplyUntransposed(vector);
        }

        private Vector MultiplyTransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
            var result = new double[NumColumns];
            // Transpose it conceptually and multiply with the vector on the right. 
            // It works because there is at most 1 non-zero per row/column
            foreach (var wholeRow in data)
            {
                int i = wholeRow.Key;
                (int j, int sign) = wholeRow.Value;
                result[j] = sign * vector[i];
            }
            return Vector.CreateFromArray(result, false);
        }

        private Vector MultiplyUntransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
            var result = new double[NumRows];
            foreach (var wholeRow in data)
            {
                int i = wholeRow.Key;
                (int j, int sign) = wholeRow.Value;
                result[i] = sign * vector[j];
            }
            return Vector.CreateFromArray(result, false);
        }
    }
}
