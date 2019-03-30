using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Sparse matrix with the non-zero entries being 1 or -1. Its main use is in domain decomposition solvers. In this context, 
    /// a <see cref="SignedBooleanMatrix"/> represents the equations that enforce continuity between the freedom degrees of 
    /// subdomains. Each row corresponds to a displacement continuity equation. If this matrix is for the whole domain, all rows 
    /// will have exactly 2 non zero entries (1 and -1). If it is only for a subdomain then some rows may be empty.
    /// Each column corresponds to a freedom degree of one of the subdomains. Columns that correspond to freedom degrees with 
    /// multiplicity = 2 will have 1 non zero entry (1 or -1). Columns that correspond to freedom degrees with multiplicity > 2, 
    /// will have at most 2 non zero entries (1 and/or -1). The other columns do not correspond to boundary dofs and will be 
    /// empty. The internal data structures that store the non-zero entries are row major.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SignedBooleanMatrix: IIndexable2D, ISparseMatrix
    {
        /// <summary>
        /// (row, (column, sign))
        /// </summary>
        private readonly Dictionary<int, Dictionary<int, int>> data;

        /// <summary>
        /// Initializes a new instance of <see cref="SignedBooleanMatrix"/> with the provided dimensions.
        /// </summary>
        /// <param name="numRows">The number of rows of the new matrix.</param>
        /// <param name="numColumns">The number of columns of the new matrix. </param>
        public SignedBooleanMatrix(int numRows, int numColumns)
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
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowIdx"/> or <paramref name="colIdx"/> violate 
        ///     the described constraints.</exception>
        public int this[int rowIdx, int colIdx]
        {
            get
            {
                if (data.TryGetValue(rowIdx, out Dictionary<int, int> colSigns))
                {
                    if (colSigns.TryGetValue(colIdx, out int sign)) return sign;
                }
                return 0;
            }
        }

        /// <summary>
        /// Sets the entry with indices (<paramref name="rowIdx"/>, <paramref name="colIdx"/>) to +1 or -1.
        /// </summary>
        /// <param name="rowIdx">The row index of the entry to set. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.</param>
        /// <param name="colIdx">The column index of the entry to set. Constraints: 
        ///     0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.</param>
        /// <param name="sign">If true, the entry (<paramref name="rowIdx"/>, <paramref name="colIdx"/>)  will be set to +1. If 
        ///     false, it will be set to -1.</param>
        public void AddEntry(int rowIdx, int colIdx, bool sign)
        {
            if (data.TryGetValue(rowIdx, out Dictionary<int, int> colSigns))
            {
                colSigns.Add(colIdx, (sign ? 1 : -1));
            }
            else
            {
                var newColSigns = new Dictionary<int, int>();
                newColSigns.Add(colIdx, (sign ? 1 : -1));
                data.Add(rowIdx, newColSigns);
            }
        }

        /// <summary>
        /// Copies the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="NumRows"/> 
        /// and length(1) = <see cref="NumColumns"/>. 
        /// </summary>
        public double[,] CopyToArray2D() => DenseStrategies.CopyToArray2D(this);

        /// <summary>
        /// Initializes a new <see cref="Matrix"/> instance by copying the entries of this <see cref="SignedBooleanMatrix"/>.
        /// </summary>
        /// <param name="transpose">If true, the new matrix will be transpose to this <see cref="SignedBooleanMatrix"/>. If 
        ///     false, thay will represent the exact same matrix (in different formats).</param>
        public Matrix CopyToFullMatrix(bool transpose)
        {
            // TODO: perhaps I should work with th col major arrays.
            if (transpose)
            {
                var dense = Matrix.CreateZero(this.NumColumns, this.NumRows);
                foreach (var wholeRow in data)
                {
                    foreach (var colValuePair in wholeRow.Value) dense[colValuePair.Key, wholeRow.Key] = colValuePair.Value;
                }
                return dense;
            }
            else
            {
                var dense = Matrix.CreateZero(this.NumRows, this.NumColumns);
                foreach (var wholeRow in data)
                {
                    foreach (var colValuePair in wholeRow.Value) dense[wholeRow.Key, colValuePair.Key] = colValuePair.Value;
                }
                return dense;
            }
        }

        /// <summary>
        /// Returns a dictionary, such that: The keys are indices of rows with at least 1 non-zero entry. The values are vectors
        /// copied from these rows.
        /// </summary>
        public Dictionary<int, Vector> CopyNonZeroRowsToVectors()
        {
            var nonZeroRows = new Dictionary<int, Vector>();
            for (int i = 0; i < NumRows; ++i)
            {
                bool exists = data.TryGetValue(i, out Dictionary<int, int> columns);
                if (exists && (columns.Count > 0)) // TODO: checking if the Dictionary<int, int> of each row is empty is probably useless. If it exists it has >= 1 entries
                {
                    var rowVector = new double[NumColumns];
                    foreach (var colValuePair in columns) rowVector[colValuePair.Key] = colValuePair.Value;
                    nonZeroRows.Add(i, Vector.CreateFromArray(rowVector));
                }
            }
            return nonZeroRows;
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.CountNonZeros"/>.
        /// </summary>
        public int CountNonZeros()
        {
            int count = 0;
            foreach (var wholeRow in data.Values) count += wholeRow.Count;
            return count;
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.EnumerateNonZeros"/>.
        /// </summary>
        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            foreach (var wholeRow in data)
            {
                foreach (var colVal in wholeRow.Value)
                {
                    yield return (wholeRow.Key, colVal.Key, colVal.Value);
                }
            }
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        /// <summary>
        /// Returns a list with the indices of the rows that have at least 1 non-zero entry (which is either 1.0 or -1.0).
        /// </summary>
        public IReadOnlyList<int> FindNonZeroRows()
        {
            // TODO: checking if the Dictionary<int, int> of each row is empty is probably useless. If it exists it has >= 1 entries

            var nonZeroRows = new List<int>();
            for (int i = 0; i < NumRows; ++i)
            {
                bool exists = data.TryGetValue(i, out Dictionary<int, int> columns);
                if (exists && (columns.Count > 0)) nonZeroRows.Add(i);
            }
            return nonZeroRows;
        }

        /// <summary>
        /// Returns a vector with the entries of the original matrix's row at index = <paramref name="rowIdx"/>.
        /// </summary>
        /// <param name="rowIdx">The index of the row to return. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="IIndexable2D.NumRows"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowIdx"/> violates the described 
        ///     constraints.</exception>
        public Vector GetRow(int rowIdx)
        {
            var rowVector = new double[NumColumns];
            bool exists = data.TryGetValue(rowIdx, out Dictionary<int, int> columns);
            if (exists)
            {
                foreach (var colValuePair in columns) rowVector[colValuePair.Key] = colValuePair.Value;
            }
            return Vector.CreateFromArray(rowVector);
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.GetSparseFormat"/>.
        /// </summary>
        public SparseFormat GetSparseFormat()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Returns true if this matrix describes the equations that enforce continuity between boundary freedom degrees that 
        /// appear in more than one subdomains. Otherwise returns false.
        /// </summary>
        public void IsContinuityEquationsMatrix()
        {
            //TODO: perhaps I should have a dedicated builder class. 
            //      Then the checks could also be more focused depending on if the matrix is global or for a subdomain
            //TODO: dedicated exception class
            foreach (var wholeRow in data)
            {
                int[] signs = wholeRow.Value.Values.ToArray();
                if (signs.Length > 2) throw new Exception(
                    $"Each row may have at most 2 entries, but row {wholeRow.Key} has {signs.Length}");
                else if ((signs.Length == 2) && (signs[0] != -signs[1])) throw new Exception(
                    $"Row {wholeRow.Key} must have two opposite signs, but they were {signs[0]} and {signs[1]}");
            }
        }

        /// <summary>
        /// Performs the matrix-vector multiplication: oper(this) * <paramref name="vector"/>.
        /// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
        /// </summary>
        /// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to the 
        ///     <see cref="IIndexable2D.NumColumns"/> of oper(this).</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
        ///     <paramref name="vector"/> is different than the <see cref="NumColumns"/> of oper(this).</exception>
        public Vector Multiply(Vector vector, bool transposeThis)
        {
            //TODO: I think that dealing with arrays will be faster than iterating the dictionaries. Another reason to separate 
            //      construction from multiplications.
            if (transposeThis) return MultiplyTransposed(vector);
            else return MultiplyUntransposed(vector);
        }

        /// <summary>
        /// Initializes a new <see cref="SignedBooleanMatrix"/> instance, that is transpose to this: result[i, j] = this[j, i]. 
        /// The entries will be explicitly copied. This method is meant for testing purposes and thus is not efficient.
        /// </summary>
        public SignedBooleanMatrix Transpose()
        {
            var transpose = new SignedBooleanMatrix(NumColumns, NumRows);
            foreach (var wholeRow in data)
            {
                foreach (var colSign in wholeRow.Value)
                {
                    transpose.AddEntry(colSign.Key, wholeRow.Key, colSign.Value == 1);
                }
            }
            return transpose;
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
                foreach (var colSign in wholeRow.Value)
                {
                    result[colSign.Key] += colSign.Value * vector[wholeRow.Key];
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
                foreach (var colSign in wholeRow.Value)
                {
                    result[wholeRow.Key] += colSign.Value * vector[colSign.Key];
                }
            }
            return Vector.CreateFromArray(result, false);
        }

        //TODO: delete this. There is now a dedicated Writer class
        //private void WriteToConsole()
        //{
        //    for (int i = 0; i < NumRows; ++i)
        //    {
        //        bool rowExists = data.TryGetValue(i, out Dictionary<int, int> colSigns);
        //        for (int j = 0; j < NumColumns; ++j)
        //        {
        //            int val = 0;
        //            if (rowExists) colSigns.TryGetValue(j, out val);
        //            Console.Write($"{val,3}");
        //        }
        //        Console.WriteLine();
        //    }
        //    Console.WriteLine();
        //}
    }
}