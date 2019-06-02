using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices.Operators
{
    /// <summary>
    /// Sparse matrix with the non-zero entries being 1 or -1. Its main use is in domain decomposition solvers. 
    /// The internal data structures that store the non-zero entries are column major.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SignedBooleanMatrixColMajor : IIndexable2D, IMappingMatrix
    {
        /// <summary>
        /// Non-zero entries: (column, (row, sign))
        /// </summary>
        private readonly Dictionary<int, Dictionary<int, int>> data;

        /// <summary>
        /// Initializes a new instance of <see cref="SignedBooleanMatrixColMajor"/> with the provided dimensions.
        /// </summary>
        /// <param name="numRows">The number of rows of the new matrix.</param>
        /// <param name="numColumns">The number of columns of the new matrix. </param>
        public SignedBooleanMatrixColMajor(int numRows, int numColumns)
        {
            this.NumRows = numRows;
            this.NumColumns = numColumns;
            this.data = new Dictionary<int, Dictionary<int, int>>();
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
                if (data.TryGetValue(colIdx, out Dictionary<int, int> rowSigns))
                {
                    if (rowSigns.TryGetValue(rowIdx, out int sign)) return sign;
                }
                return 0;
            }
        }

        /// <summary>
        /// Sets the entry with indices (<paramref name="rowIdx"/>, <paramref name="colIdx"/>) to +1 or -1.
        /// </summary>
        /// <param name="rowIdx">
        /// The row index of the entry to set. Constraints: 
        /// 0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.
        /// </param>
        /// <param name="colIdx">
        /// The column index of the entry to set. Constraints: 
        /// 0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.
        /// </param>
        /// <param name="sign">
        /// If true, the entry (<paramref name="rowIdx"/>, <paramref name="colIdx"/>)  will be set to +1. 
        /// If false, it will be set to -1.
        /// </param>
        public void AddEntry(int rowIdx, int colIdx, bool sign)
        {
            if (data.TryGetValue(colIdx, out Dictionary<int, int> rowSigns))
            {
                rowSigns.Add(rowIdx, (sign ? 1 : -1));
            }
            else
            {
                var newRowSigns = new Dictionary<int, int>();
                newRowSigns.Add(rowIdx, (sign ? 1 : -1));
                data.Add(colIdx, newRowSigns);
            }
        }

        /// <summary>
        /// Initializes a new <see cref="Matrix"/> instance by copying the entries of this <see cref="SignedBooleanMatrixRowMajor"/>.
        /// </summary>
        /// <param name="transpose">
        /// If true, the new matrix will be transpose to this <see cref="SignedBooleanMatrixRowMajor"/>. If false, they will 
        /// represent the exact same matrix (in different formats).
        /// </param>
        public Matrix CopyToFullMatrix(bool transpose)
        {
            // TODO: perhaps I should work with th col major arrays.
            if (transpose)
            {
                var dense = Matrix.CreateZero(this.NumColumns, this.NumRows);
                foreach (var wholeCol in data)
                {
                    foreach (var rowValuePair in wholeCol.Value) dense[wholeCol.Key, rowValuePair.Key] = rowValuePair.Value;
                }
                return dense;
            }
            else
            {
                var dense = Matrix.CreateZero(this.NumRows, this.NumColumns);
                foreach (var wholeCol in data)
                {
                    foreach (var rowValuePair in wholeCol.Value) dense[rowValuePair.Key, wholeCol.Key] = rowValuePair.Value;  
                }
                return dense;
            }
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        public SignedBooleanMatrixColMajor GetColumns(int[] colsToKeep, bool deepCopy)
        {
            //TODO: This would be more efficient if the matrix was column major
            var clone = new SignedBooleanMatrixColMajor(this.NumRows, colsToKeep.Length);
            if (deepCopy)
            {
                for (int j = 0; j < colsToKeep.Length; ++j)
                {
                    var cloneColumn = new Dictionary<int, int>();
                    foreach (var rowSign in this.data[colsToKeep[j]])
                    {
                        cloneColumn[rowSign.Key] = rowSign.Value;
                    }
                    clone.data[j] = cloneColumn;
                }
            }
            else
            {
                for (int j = 0; j < colsToKeep.Length; ++j)
                {
                    clone.data[j] = this.data[colsToKeep[j]];
                }
            }
            return clone;
        }

        public Vector Multiply(Vector vector, bool transposeThis = false)
        {
            //TODO: I think that dealing with arrays will be faster than iterating the dictionaries. Another reason to separate 
            //      construction from multiplications.
            if (transposeThis) return MultiplyTransposed(vector);
            else return MultiplyUntransposed(vector);
        }

        public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
        {
            if (transposeThis) return MultiplyRightTransposed(other);
            else return MultiplyRightUntransposed(other);
        }

        private Vector MultiplyTransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
            var result = new double[NumColumns];
            foreach (var wholeCol in data)
            {
                double sum = 0.0;
                foreach (var rowSign in wholeCol.Value)
                {
                    sum += rowSign.Value * vector[rowSign.Key];
                }
                result[wholeCol.Key] = sum;
            }
            return Vector.CreateFromArray(result, false);
        }

        private Vector MultiplyUntransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
            var result = new double[NumRows];
            // Transpose it conceptually and multiply with the vector on the right. 
            foreach (var wholeCol in data)
            {
                foreach (var rowSign in wholeCol.Value)
                {
                    result[rowSign.Key] += rowSign.Value * vector[wholeCol.Key];
                }
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
                int offset = j * this.NumColumns;
                foreach (var wholeCol in data)
                {
                    double sum = 0.0;
                    foreach (var rowSign in wholeCol.Value)
                    {
                        sum += rowSign.Value * other[rowSign.Key, j];
                    }
                    result[offset + wholeCol.Key] = sum;
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
                foreach (var wholeCol in data)
                {
                    foreach (var rowSign in wholeCol.Value)
                    {
                        result[offset + rowSign.Key] += rowSign.Value * other[wholeCol.Key, j];
                    }
                }
            }
            return Matrix.CreateFromArray(result, this.NumRows, other.NumColumns, false);
        }
    }
}
