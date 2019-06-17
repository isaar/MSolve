using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
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
    public class UnsignedBooleanMatrix : IIndexable2D, IMappingMatrix
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

        /// <summary>
        /// WARNING: this only works if this matrix is a mapping matrix L used in FETI solvers, meaning :
        /// 1) There are more columns than rows.
        /// 2) Each row of this matrix must have exactly one 1 and all other 0,
        /// 3) Each column must have at most one 1. It is possible that a column is completely 0.
        /// </summary>
        //TODO: This should be the actual way this matrix is stored. The dictionaries should be for building it only.
        public int[] GetRowsToColumnsMap()
        {
            var map = new int[NumRows];
            for (int i = 0; i < NumRows; ++i)
            {
                Debug.Assert(data[i].Count == 1);
                map[i] = data[i].First();
            }
            return map;
        }

        public Vector Multiply(Vector vector, bool transposeThis = false)
        {
            if (transposeThis) return MultiplyTransposed(vector);
            else return MultiplyUntransposed(vector);
        }

        public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
        {
            if (transposeThis) return MultiplyRightTransposed(other);
            else return MultiplyRightUntransposed(other);
        }

        /// <summary>
        /// WARNING: this only works if this matrix is a mapping matrix L used in FETI solvers, meaning :
        /// 1) There are more columns than rows.
        /// 2) Each row of this matrix must have exactly one 1 and all other 0,
        /// 3) Each column must have at most one 1. It is possible that a column is completely 0.
        /// Essentially this method performs global matrix assembly: globalK = L^T * localK * L.
        /// </summary>
        /// <param name="other"></param>
        public Matrix ThisTransposeTimesOtherTimesThis(Matrix other) //TODO: this should be implemented for symmetric matrices
        {
            //TODO: Move this class to project Solvers where the following assumptions are always satisfied.
            //TODO: Otherwise, rename this method to specify for which instances it works correctly.

            // Rows of this matrix correspond to rows of matrix "other" (local). Columns of this matrix correspond to columns
            // of matrix "result" (global).
            Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumRows);
            var result = Matrix.CreateZero(this.NumColumns, this.NumColumns);
            for (int otherCol = 0; otherCol < other.NumColumns; ++otherCol)
            {
                Debug.Assert(data[otherCol].Count == 1);
                int resultCol = data[otherCol].First();
                for (int otherRow = 0; otherRow < other.NumRows; ++otherRow)
                {
                    Debug.Assert(data[otherRow].Count == 1);
                    int resultRow = data[otherRow].First();
                    result[resultRow, resultCol] = other[otherRow, otherCol];
                }
            }
            return result;
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

        private Matrix MultiplyRightTransposed(Matrix other)
        {
            //TODO: I think that it will pay off to transpose an all integer CSR matrix and store both. Especially in the case 
            //     of subdomain boolean matrices, that little extra memory should not be of concern.
            Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
            var result = new double[this.NumColumns * other.NumRows];
            for (int j = 0; j < other.NumColumns; ++j)
            {
                int offset = j * this.NumRows;
                // Transpose it conceptually and multiply with the vector on the right. 
                foreach (var wholeRow in data)
                {
                    foreach (int col in wholeRow.Value)
                    {
                        result[offset + col] += other[wholeRow.Key, j];
                    }
                }
            }
            return Matrix.CreateFromArray(result, this.NumColumns, other.NumColumns, false);
        }

        private Matrix MultiplyRightUntransposed(Matrix other)
        {
            Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
            var result = new double[this.NumRows * other.NumColumns];
            for (int j = 0; j < other.NumColumns; ++j)
            {
                int offset = j * this.NumRows;
                foreach (var wholeRow in data)
                {
                    double sum = 0.0;
                    foreach (int col in wholeRow.Value)
                    {
                        sum += other[col, j];
                    }
                    result[offset + wholeRow.Key] = sum;
                }
            }
            return Matrix.CreateFromArray(result, this.NumRows, other.NumColumns, false);
        }
    }
}
