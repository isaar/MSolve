using System;
using System.Collections.Generic;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

// TODO: try to make general versions of row major and col major multiplication. The lhs matrix/vector will be supplied by the 
// caller, depending on if it is a matrix, a transposed matrix, a vector, etc. Compare its performance with the verbose code and 
// also try inlining. C preprocessor macros would be actually useful here. Otherwise, move all that boilerplate code to a 
// CSRStrategies static class.
// TODO: In matrix-matrix/vector multiplications: perhaps I should work with a column major array directly instead of an output   
// Matrix and an array instead of an output Vector.
// TODO: perhaps optimizations if (other is Matrix) are needed, to directly index into its raw col major array.
// The access paterns are always the same
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Sparse matrix stored in Compressed Sparse Columns format (3-array version). The CSR format is optimized for matrix-vector 
    /// and matrix-matrix multiplications, where the CSC matrix is on the left transposed or on the right untransposed. The other
    /// multiplicationss are more efficient using <see cref="CsrMatrix"/>. To build a <see cref="CscMatrix"/> conveniently, 
    /// use <see cref="Builders.DokColMajor"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CscMatrix: IMatrix, ISparseMatrix //TODO: Use MKL with descriptors
    {
        private readonly double[] values;
        private readonly int[] rowIndices;
        private readonly int[] colOffsets;

        private CscMatrix(int numRows, int numCols, double[] values, int[] rowIndices, int[] colOffsets)
        {
            this.values = values;
            this.rowIndices = rowIndices;
            this.colOffsets = colOffsets;
            this.NumRows = numRows;
            this.NumColumns = numCols;
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
        public double this[int rowIdx, int colIdx]
        {
            get
            {
                int index = FindIndexOf(rowIdx, colIdx);
                if (index == -1) return 0.0;
                else return values[index];
            }
        }

        /// <summary>
        /// Initializes a new <see cref="CscMatrix"/> with the specified dimensions and the provided arrays 
        /// (<paramref name="values"/>, <paramref name="rowIndices"/> and <paramref name="colOffsets"/>) as its internal data.
        /// </summary>
        /// <param name="numRows">The number of rows of the new matrix.</param>
        /// <param name="numCols">The number of columns of the new matrix.</param>
        /// <param name="values">Array that contains the non-zero entries. It must have the same length 
        ///     as <paramref name="rowIndices"/>. The non-zero entries of each column must appear consecutively in 
        ///     <paramref name="values"/>. They can also be sorted in increasing order of their row indices, which speeds up
        ///     subsequent operations.</param>
        /// <param name="rowIndices">Array that contains the row indices of the non-zero entries. It must have the same 
        ///     length as <paramref name="values"/>. There is an 1 to 1 matching between these two arrays: 
        ///     <paramref name="rowIndices"/>[i] is the row index of the entry <paramref name="values"/>[i]. Also:
        ///     0 &lt;= <paramref name="rowIndices"/>[i] &lt; <paramref name="numRows"/>.</param>
        /// <param name="colOffsets">Array that contains the index of the first entry of each column into the arrays 
        ///     <paramref name="values"/> and <paramref name="rowIndices"/>. Its length is <paramref name="numRows"/> + 1. The 
        ///     last entry is the number of non-zero entries, which must be equal to the length of <paramref name="values"/> 
        ///     and <paramref name="rowIndices"/>.</param>
        /// <param name="checkInput">If true, the provided arrays will be checked to make sure they are valid CSC arrays, which 
        ///     is safer. If false, no such check will take place, which is faster.</param>
        public static CscMatrix CreateFromArrays(int numRows, int numCols, double[] values, int[] rowIndices, int[] colOffsets,
            bool checkInput)
        {
            if (checkInput)
            {
                if (colOffsets.Length != numCols + 1)
                {
                    throw new ArgumentException("The length of the CSC column offsets array must be equal to the number of"
                        + " rows + 1, but was " + colOffsets.Length);
                }
                if (values.Length != rowIndices.Length)
                {
                    throw new ArgumentException("The length of the CSC values and row indices arrays must be equal (and equal"
                        + $" to the number of non zero entries), but were {values.Length} and {rowIndices.Length} respectively");
                }
                if (colOffsets[0] != 0)
                {
                    throw new ArgumentException("The first entry of the CSC column offsets array must be 0, but was " 
                        + colOffsets[0]);
                }
                if (colOffsets[colOffsets.Length - 1] != values.Length)
                {
                    throw new ArgumentException("The last entry of the CSC column offsets array must be equal to the number of"
                        + " non zero entries, but was " + colOffsets[colOffsets.Length - 1]);
                }
            }
            return new CscMatrix(numRows, numCols, values, rowIndices, colOffsets);
        }

        #region operators (use extension operators when they become available)
        /// <summary>
        /// Performs the matrix-vector multiplication: result = <paramref name="vectorLeft"/> * <paramref name="matrixRight"/>.
        /// If <paramref name="matrixRight"/> is m1-by-n1 and <paramref name="vectorLeft"/> has length = n2, then m1 must be 
        /// equal to n2. The result will be a vector with length = n1, written to a new <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="vectorLeft">The <see cref="Vector"/> operand on the left. It can be considered as a row vector.</param>
        /// <param name="matrixRight">The <see cref="CscMatrix"/> operand on the right.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrixRight"/>.<see cref="NumRows"/> is 
        ///     different than <paramref name="vectorLeft"/>.<see cref="Vector.Length"/>.</exception>
        public static Vector operator *(Vector vectorLeft, CscMatrix matrixRight)
            => matrixRight.MultiplyRight(vectorLeft, true);
        #endregion

        /// <summary>
        /// See <see cref="IMatrixView.Axpy(IMatrixView, double)"/>.
        /// </summary>
        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CscMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    Array.Copy(this.values, resultValues, values.Length);
                    CBlas.Daxpy(values.Length, otherCoefficient, ref otherCSC.values[0], 1, ref resultValues[0], 1);
                    return new CscMatrix(NumRows, NumColumns, resultValues, this.rowIndices, this.colOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, 1.0, otherMatrix, otherCoefficient);
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix is written to a new <see cref="CscMatrix"/> and then returned.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
        ///     <see cref="CscMatrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
        public CscMatrix Axpy(CscMatrix otherMatrix, double otherCoefficient)
        {
            // Conceptually it is not wrong to so this, even if the indexers are different, but how would I implement it.
            if (!HasSameIndexer(otherMatrix))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] resultValues = new double[values.Length];
            Array.Copy(this.values, resultValues, values.Length);
            CBlas.Daxpy(values.Length, otherCoefficient, ref otherMatrix.values[0], 1, ref resultValues[0], 1);
            // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
            return new CscMatrix(NumRows, NumColumns, resultValues, this.rowIndices, this.colOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrix.AxpyIntoThis(IMatrixView, double)"/>.
        /// </summary>
        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CscMatrix casted) AxpyIntoThis(casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                 "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="CscMatrix"/> instance.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same indexing arrays as this <see cref="CscMatrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="SparsityPatternModifiedException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     indexing arrays than this instance.</exception>
        public void AxpyIntoThis(CscMatrix otherMatrix, double otherCoefficient)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (!HasSameIndexer(otherMatrix))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            CBlas.Daxpy(values.Length, otherCoefficient, ref otherMatrix.values[0], 1, ref this.values[0], 1);
        }

        /// <summary>
        /// Initializes a new <see cref="Matrix"/> instance by copying the entries of this <see cref="CsrMatrix"/>. 
        /// </summary>
        public Matrix CopyToFullMatrix()
        {
            Matrix fullMatrix = Matrix.CreateZero(this.NumRows, this.NumColumns);
            for (int j = 0; j < this.NumColumns; ++j) //Row major order
            {
                int colStart = colOffsets[j]; //inclusive
                int colEnd = colOffsets[j + 1]; //exclusive
                for (int k = colStart; k < colEnd; ++k)
                {
                    fullMatrix[rowIndices[k], j] = values[k];
                }
            }
            return fullMatrix;
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.CountNonZeros"/>
        /// </summary>
        public int CountNonZeros() => values.Length;

        /// <summary>
        /// See <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public IMatrixView DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is CscMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherCSC.values[i]);
                    }
                    return new CscMatrix(NumRows, NumColumns, resultValues, rowIndices, colOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.DoEntrywise(this, other, binaryOperation);
        }

        /// <summary>
        /// See <see cref="IMatrix.DoEntrywiseIntoThis(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is CscMatrix casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// this[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="matrix"/>[i, j]). 
        /// The resulting matrix overwrites the entries of this <see cref="CscMatrix"/> instance.
        /// </summary>
        /// <param name="matrix">A matrix with the same indexing arrays as this <see cref="CscMatrix"/> instance.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="SparsityPatternModifiedException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     indexing arrays than this instance.</exception>
        public void DoEntrywiseIntoThis(CscMatrix other, Func<double, double, double> binaryOperation)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (!HasSameIndexer(other))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            for (int i = 0; i < values.Length; ++i) this.values[i] = binaryOperation(this.values[i], other.values[i]);
        }

        /// <summary>
        /// See <see cref="IMatrixView.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        IMatrixView IMatrixView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            // Only apply the operation on non zero entries
            double[] newValues = new double[values.Length];
            for (int i = 0; i < values.Length; ++i) newValues[i] = unaryOperation(values[i]);

            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Copy the index arrays. TODO: See if we can use the same index arrays (e.g. if this class does not change them (it shouldn't))
                int[] rowIndicesCopy = new int[rowIndices.Length];
                Array.Copy(rowIndices, rowIndicesCopy, rowIndices.Length);
                int[] colOffsetsCopy = new int[colOffsets.Length];
                Array.Copy(colOffsets, colOffsetsCopy, colOffsets.Length);
                return new CscMatrix(NumRows, NumColumns, newValues, rowIndicesCopy, colOffsetsCopy);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                return new CscMatrix(NumRows, NumColumns, newValues, rowIndices, colOffsets).CopyToFullMatrix();
            }
        }

        /// <summary>
        /// See <see cref="IMatrix.DoToAllEntriesIntoThis(Func{double, double})"/>.
        /// </summary>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0))
            {
                for (int i = 0; i < values.Length; ++i) values[i] = unaryOperation(values[i]);
            }
            else
            {
                throw new SparsityPatternModifiedException("This operation will change the sparsity pattern");
            }
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.EnumerateNonZeros"/>.
        /// </summary>
        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            for (int j = 0; j < NumColumns; ++j)
            {
                int colStart = colOffsets[j]; //inclusive
                int colEnd = colOffsets[j + 1]; //exclusive
                for (int k = colStart; k < colEnd; ++k)
                {
                    yield return (rowIndices[k], j, values[k]);
                }
            }
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
            var comparer = new ValueComparer(1e-13);
            for (int j = 0; j < NumColumns; ++j)
            {
                int colStart = colOffsets[j]; // Inclusive
                int colEnd = colOffsets[j + 1]; // Exclusive
                int previousRow = 0;
                for (int k = colStart; k < colEnd; ++k)
                {
                    int row = rowIndices[k];
                    for (int i = previousRow; i < row; ++i) // Zero entries between the stored ones
                    {
                        if (!comparer.AreEqual(0.0, other[i, j])) return false;
                    }
                    if (!comparer.AreEqual(values[k], other[row, j])) return false; // Non zero entry
                    previousRow = row + 1;
                }
            }
            return true; // At this point all entries have been checked and are equal
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.GetSparseFormat"/>.
        /// </summary>
        public SparseFormat GetSparseFormat()
        {
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Row indices", rowIndices);
            format.RawIndexArrays.Add("Column offsets", colOffsets);
            return format;
        }

        /// <summary>
        /// See <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>.
        /// </summary>
        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CscMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    Array.Copy(this.values, resultValues, values.Length);
                    CBlas.Daxpby(values.Length, otherCoefficient, ref otherCSC.values[0], 1, thisCoefficient, ref this.values[0], 1);
                    return new CscMatrix(NumRows, NumColumns, resultValues, this.rowIndices, this.colOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, thisCoefficient, otherMatrix, otherCoefficient);
        }

        /// <summary>
        /// See <see cref="IMatrix.LinearCombinationIntoThis(double, IMatrixView, double)"/>.
        /// </summary>
        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CscMatrix casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// this[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
        ///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="CscMatrix"/> instance.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Matrix"/>.</param>
        /// <param name="otherMatrix">A matrix with the same indexing arrays as this <see cref="CscMatrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="SparsityPatternModifiedException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     indexing arrays than this instance.</exception>
        public void LinearCombinationIntoThis(double thisCoefficient, CscMatrix otherMatrix, double otherCoefficient)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (!HasSameIndexer(otherMatrix))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            CBlas.Daxpby(values.Length, otherCoefficient, ref otherMatrix.values[0], 1, thisCoefficient, ref this.values[0], 1);
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyLeft(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            
            if (transposeThis) // Compute one output row at a time. TODO: not efficent for the col major output matrix.
            {
                if (transposeOther)
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumColumns);
                    var result = Matrix.CreateZero(other.NumColumns, this.NumRows);
                    for (int r = 0; r < result.NumRows; ++r) // Compute one output row at a time.
                    {
                        // x * A^T = linear combination of rows of A^T = columns of A, with the entries of x as coefficients, 
                        // where x is row r of transpose(other matrix)
                        for (int j = 0; j < this.NumColumns; ++j)
                        {
                            double scalar = other[j, r];
                            int cscColStart = colOffsets[j]; //inclusive
                            int cscColEnd = colOffsets[j + 1]; //exclusive
                            for (int k = cscColStart; k < cscColEnd; ++k)
                            {
                                result[r, rowIndices[k]] += scalar * values[k]; // sum(other[j,r] * transpose(csc).row[j])
                            }
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumColumns);
                    var result = Matrix.CreateZero(other.NumRows, this.NumRows);
                    for (int r = 0; r < result.NumRows; ++r) // Compute one output row at a time.
                    {
                        // x * A^T = linear combination of rows of A^T = columns of A, with the entries of x as coefficients, 
                        // where x is row r of the other matrix
                        for (int j = 0; j < this.NumColumns; ++j)
                        {
                            double scalar = other[r, j];
                            int cscColStart = colOffsets[j]; //inclusive
                            int cscColEnd = colOffsets[j + 1]; //exclusive
                            for (int k = cscColStart; k < cscColEnd; ++k) // sum(other[r,j] * transpose(csc).row[j])
                            {
                                result[r, rowIndices[k]] += scalar * values[k];
                            }
                        }
                    }
                    return result;
                }
            }
            else
            {
                if (transposeOther)
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumRows);
                    var result = Matrix.CreateZero(other.NumColumns, this.NumColumns);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time. TODO: perhaps 1 row at a time.
                    {
                        int cscColStart = colOffsets[c]; //inclusive
                        int cscColEnd = colOffsets[c + 1]; //exclusive
                        for (int i = 0; i < other.NumColumns; ++i)
                        {
                            double dot = 0.0;
                            for (int k = cscColStart; k < cscColEnd; ++k) // other.col[i] * csc.col[c]
                            {
                                dot += values[k] * other[rowIndices[k], i]; 
                            }
                            result[i, c] = dot;
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumRows);
                    var result = Matrix.CreateZero(other.NumRows, this.NumColumns);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time. TODO: perhaps 1 row at a time.
                    {
                        int cscColStart = colOffsets[c]; //inclusive
                        int cscColEnd = colOffsets[c + 1]; //exclusive
                        for (int i = 0; i < other.NumRows; ++i)
                        {
                            double dot = 0.0;
                            for (int k = cscColStart; k < cscColEnd; ++k) // other.row[i] * csc.col[c]
                            {
                                dot += values[k] * other[i, rowIndices[k]];
                            }
                            result[i, c] = dot;
                        }
                    }
                    return result;
                }
            }
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            // TODO: Throwing exceptions when csc is on the left seems attractive.
            if (transposeThis)
            {
                if (transposeOther)
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumColumns);
                    var result = Matrix.CreateZero(this.NumColumns, other.NumRows);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
                    {
                        for (int j = 0; j < this.NumColumns; ++j)
                        {
                            int cscColStart = colOffsets[j]; //inclusive
                            int cscColEnd = colOffsets[j + 1]; //exclusive
                            double dot = 0.0;
                            for (int k = cscColStart; k < cscColEnd; ++k) // traspose(csc).row[i] * other.col[c]
                            {
                                dot += values[k] * other[c, rowIndices[k]];
                            }
                            result[j, c] = dot;
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
                    var result = Matrix.CreateZero(this.NumColumns, other.NumColumns);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
                    {
                        for (int j = 0; j < this.NumColumns; ++j)
                        {
                            int cscColStart = colOffsets[j]; //inclusive
                            int cscColEnd = colOffsets[j + 1]; //exclusive
                            double dot = 0.0;
                            for (int k = cscColStart; k < cscColEnd; ++k) // traspose(csc).row[i] * other.col[c]
                            {
                                dot += values[k] * other[rowIndices[k], c];
                            }
                            result[j, c] = dot;
                        }
                    }
                    return result;
                }
            }
            else
            {
                if (transposeOther)
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumColumns);
                    var result = Matrix.CreateZero(this.NumRows, other.NumRows);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
                    {
                        // A * x = linear combination of columns of A, with the entries of x as coefficients, 
                        // where x is column c of transpose(other matrix)
                        for (int j = 0; j < NumColumns; ++j)
                        {
                            double scalar = other[c, j];
                            int cscColStart = colOffsets[j]; //inclusive
                            int cscColEnd = colOffsets[j + 1]; //exclusive
                            for (int k = cscColStart; k < cscColEnd; ++k) // sum(other[c,j] * csc.col[j]))
                            {
                                result[rowIndices[k], c] += scalar * values[k];
                            }
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
                    var result = Matrix.CreateZero(this.NumRows, other.NumColumns);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
                    {
                        // A * x = linear combination of columns of A, with the entries of x as coefficients, 
                        // where x is column c of the other matrix
                        for (int j = 0; j < NumColumns; ++j)
                        {
                            double scalar = other[j, c];
                            int cscColStart = colOffsets[j]; //inclusive
                            int cscColEnd = colOffsets[j + 1]; //exclusive
                            for (int k = cscColStart; k < cscColEnd; ++k) // sum(other[j,c] * csc.col[j]))
                            {
                                result[rowIndices[k], c] += scalar * values[k];
                            }
                        }
                    }
                    return result;
                }
            }
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IVectorView, bool)"/>.
        /// </summary>
        public Vector MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            if (vector is Vector) return MultiplyRight((Vector)vector, transposeThis);
            else throw new NotImplementedException();
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
        public Vector MultiplyRight(Vector vector, bool transposeThis = false)
        {
            if (transposeThis)
            {
                // A^T * x = sum{row of A^T * x} = sum{col of A * x}
                Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
                double[] result = new double[NumColumns];
                for (int j = 0; j < NumColumns; ++j)
                {
                    double dot = 0.0;
                    int colStart = colOffsets[j]; //inclusive
                    int colEnd = colOffsets[j + 1]; //exclusive
                    for (int k = colStart; k < colEnd; ++k)
                    {
                        dot += values[k] * vector[rowIndices[k]];
                    }
                    result[j] = dot;
                }
                return Vector.CreateFromArray(result, false);
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
                // A * x = linear combination of columns of A, with the entries of x as coefficients
                double[] result = new double[NumRows];
                for (int j = 0; j < NumColumns; ++j)
                {
                    double scalar = vector[j];
                    int colStart = colOffsets[j]; //inclusive
                    int colEnd = colOffsets[j + 1]; //exclusive
                    for (int k = colStart; k < colEnd; ++k)
                    {
                        result[rowIndices[k]] += scalar * values[k];
                    }
                }
                return Vector.CreateFromArray(result, false);
            }
        }

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            int nnz = values.Length;
            for (int i = 0; i < nnz; ++i) aggregator = processEntry(values[i], aggregator);
            aggregator = processZeros(NumRows * NumColumns - nnz, aggregator);
            return finalize(aggregator);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Scale(double)"/>.
        /// </summary>
        IMatrixView IMatrixView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// Performs the following operation for the non-zero entries (i, j), such that 0 &lt;= i &lt; <see cref="NumRows"/>, 
        /// 0 &lt;= j &lt; <see cref="NumColumns"/>: result[i, j] = <paramref name="scalar"/> * this[i, j].
        /// The resulting matrix is written to a new <see cref="CscMatrix"/> and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
        public CscMatrix Scale(double scalar)
        {
            int nnz = this.values.Length;
            double[] resultValues = new double[nnz]; 
            Array.Copy(this.values, resultValues, nnz); //TODO: perhaps I should also copy the indexers
            CBlas.Dscal(nnz, scalar, ref resultValues[0], 1);
            return new CscMatrix(this.NumRows, this.NumColumns, resultValues, this.rowIndices, this.colOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrix.ScaleIntoThis(double)"/>.
        /// </summary>
        public void ScaleIntoThis(double scalar) => CBlas.Dscal(values.Length, scalar, ref values[0], 1);

        /// <summary>
        /// See <see cref="IMatrix.SetEntryRespectingPattern(int, int, double)"/>.
        /// </summary>
        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            int index = FindIndexOf(rowIdx, colIdx);
            if (index == -1) throw new SparsityPatternModifiedException($"Cannot write to zero entry ({rowIdx}, {colIdx}).");
            else values[index] = value;
        }

        /// <summary>
        /// See <see cref="IMatrixView.Transpose"/>.
        /// </summary>
        public IMatrixView Transpose() => TransposeToCSR(true);

        /// <summary>
        /// Creates a new <see cref="CscMatrix"/> instance, that is transpose to this: result[i, j] = this[j, i].
        /// </summary>
        public CscMatrix TransposeToCSC()
        {
            // Use C# port of the scipy method.
            // TODO: Perhaps it could be done faster by making extra assumptions. Otherwise use MKL
            int nnz = this.values.Length;
            var csrValues = new double[nnz];
            var csrColIndices = new int[nnz];
            var csrRowOffsets = new int[NumRows + 1];

            SparseArrays.CsrToCsc(NumColumns, NumRows, this.colOffsets, this.rowIndices, this.values,
                csrRowOffsets, csrColIndices, csrValues);

            return new CscMatrix(NumColumns, NumRows, csrValues, csrColIndices, csrRowOffsets);
        }

        /// <summary>
        /// Creates a new <see cref="CsrMatrix"/> instance, that is transpose to this: result[i, j] = this[j, i]. The 
        /// internal arrays can be copied or shared with this <see cref="CscMatrix"/> instance.
        /// </summary>
        /// <param name="copyInternalArray">If true, the internal arrays that store the entries of this 
        ///     <see cref="CscMatrix"/> instance will be copied and the new <see cref="CsrMatrix"/> instance 
        ///     instance will have references to the copies, which is safer. If false, both the new matrix and this one will have  
        ///     references to the same internal arrays, which is faster.</param>
        public CsrMatrix TransposeToCSR(bool copyInternalArrays)
        {
            if (copyInternalArrays)
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                int[] rowIndicesCopy = new int[rowIndices.Length];
                Array.Copy(rowIndices, rowIndicesCopy, rowIndices.Length);
                int[] colOffsetsCopy = new int[colOffsets.Length];
                Array.Copy(colOffsets, colOffsetsCopy, colOffsets.Length);
                return CsrMatrix.CreateFromArrays(NumColumns, NumRows, valuesCopy, rowIndicesCopy, colOffsetsCopy, false);
            }
            else return CsrMatrix.CreateFromArrays(NumColumns, NumRows, values, rowIndices, colOffsets, false);
        }

        /// <summary>
        /// Return the index into values and rowIndices arrays, if the (rowIdx, colIdx) entry is within the pattern. 
        /// Otherwise returns -1.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        /// <returns></returns>
        private int FindIndexOf(int rowIdx, int colIdx)
        {
            Preconditions.CheckIndices(this, rowIdx, colIdx); //TODO: check indices?
            int colStart = colOffsets[colIdx]; //inclusive
            int colEnd = colOffsets[colIdx + 1]; //exclusive
            for (int k = colStart; k < colEnd; ++k)
            {
                if (rowIndices[k] == rowIdx) return k;
            }
            return -1;
        }

        private bool HasSameIndexer(CscMatrix other)
        {
            return (this.rowIndices == other.rowIndices) && (this.colOffsets == other.colOffsets);
        }
    }
}
