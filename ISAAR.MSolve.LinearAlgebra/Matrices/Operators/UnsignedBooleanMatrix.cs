using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices.Operators
{
    /// <summary>
    /// Sparse matrix with the non-zero entries being 1. Its main use is in domain decomposition solvers. 
    /// The internal data structures that store the non-zero entries are row major.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class UnsignedBooleanMatrix : IMappingMatrix
    {
        /// <summary>
        /// Non-zero entries: (row, columns)
        /// </summary>
        private readonly Dictionary<int, HashSet<int>> data;

        /// <summary>
        /// Initializes a new instance of <see cref="UnsignedBooleanMatrix"/> with the provided dimensions.
        /// </summary>
        /// <param name="numRows">The number of rows of the new matrix.</param>
        /// <param name="numColumns">The number of columns of the new matrix. </param>
        public UnsignedBooleanMatrix(int numRows, int numColumns)
        {
            this.NumRows = numRows;
            this.NumColumns = numColumns;
            this.data = new Dictionary<int, HashSet<int>>();
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>.
        /// </summary>
        /// <remarks>
        /// The entries can be 0.0, 1.0 or -1.0
        /// </remarks>
        double IIndexable2D.this[int rowIdx, int colIdx] => this[rowIdx, colIdx];

        /// <summary>
        /// The entry with row index = rowIdx and column index = colIdx. 
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.</param>
        /// <param name="colIdx">The column index: 0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.</param>
        /// <exception cref="IndexOutOfRangeException">
        /// Thrown if <paramref name="rowIdx"/> or <paramref name="colIdx"/> violate the described constraints.
        /// </exception>
        public int this[int rowIdx, int colIdx]
        {
            get
            {
                if (data.TryGetValue(rowIdx, out HashSet<int> colIndices))
                {
                    if (colIndices.Contains(colIdx)) return 1;
                }
                return 0;
            }
        }

        /// <summary>
        /// Sets the entry with indices (<paramref name="rowIdx"/>, <paramref name="colIdx"/>) to 1.
        /// </summary>
        /// <param name="rowIdx">
        /// The row index of the entry to set. Constraints: 
        /// 0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.
        /// </param>
        /// <param name="colIdx">
        /// The column index of the entry to set. Constraints: 
        /// 0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.
        /// </param>
        public void AddEntry(int rowIdx, int colIdx)
        {
            if (data.TryGetValue(rowIdx, out HashSet<int> colIndices))
            {
                colIndices.Add(colIdx);
            }
            else
            {
                colIndices = new HashSet<int>();
                data[rowIdx] = colIndices;
            }
            colIndices.Add(colIdx);
        }

        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        public Vector Multiply(Vector vector, bool transposeThis = true)
        {
            if (transposeThis) return MultiplyTransposed(vector);
            else return MultiplyUntransposed(vector);
        }

        private Vector MultiplyTransposed(Vector vector)
        {
            //TODO: I think that it will pay off to transpose an all integer CSR matrix and store both. Especially in the case 
            //     of subdomain boolean matrices, that little extra memory should not be of concern.
            Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
            var result = new double[NumColumns];
            // Transpose it conceptually and multiply with the vector on the right. 
            foreach (var wholeRow in data)
            {
                foreach (int col in wholeRow.Value)
                {
                    result[col] += vector[wholeRow.Key];
                }
            }
            return Vector.CreateFromArray(result, false);
        }

        private Vector MultiplyUntransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
            var result = new double[NumRows];
            foreach (var wholeRow in data)
            {
                double sum = 0.0;
                foreach (int col in wholeRow.Value)
                {
                    sum += vector[col];
                }
                result[wholeRow.Key] = sum;
            }
            return Vector.CreateFromArray(result, false);
        }
    }
}
