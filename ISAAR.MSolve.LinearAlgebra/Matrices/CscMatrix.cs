using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using static ISAAR.MSolve.LinearAlgebra.LibrarySettings;

//TODO: try to make general versions of row major and col major multiplication. The lhs matrix/vector will be supplied by the 
// caller, depending on if it is a matrix, a transposed matrix, a vector, etc. Compare its performance with the verbose code and 
// also try inlining. C preprocessor macros would be actually useful here. Otherwise, move all that boilerplate code to a 
// CSRStrategies static class.
//TODO: In matrix-matrix/vector multiplications: perhaps I should work with a column major array directly instead of an output   
//      Matrix and an array instead of an output Vector.
//TODO: perhaps optimizations if (other is Matrix) are needed, to directly index into its raw col major array.
//      The access paterns are always the same.
//TODO: The implementations of this class should call transposed operations on a backing CSR matrix.
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Sparse matrix stored in Compressed Sparse Columns format (3-array version). The CSR format is optimized for matrix-vector 
    /// and matrix-matrix multiplications, where the CSC matrix is on the left transposed or on the right untransposed. The other
    /// multiplicationss are more efficient using <see cref="CsrMatrix"/>. To build a <see cref="CscMatrix"/> conveniently, 
    /// use <see cref="Builders.DokColMajor"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CscMatrix: IMatrix, ISparseMatrix
    {
        private const int zeroEntryOffset = -1;

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
        /// The number of non zero entries of the matrix.
        /// </summary>
        public int NumNonZeros => colOffsets[colOffsets.Length - 1];

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
                int entryOffset = FindOffsetOf(rowIdx, colIdx);
                if (entryOffset == zeroEntryOffset) return 0.0;
                else return values[entryOffset];
            }
        }

        /// <summary>
        /// The internal array that stores the non-zero entries of the matrix. The non-zero entries of each column are 
        /// consecutive. Its length is equal to the number of non-zero entries. 
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal double[] RawValues => values;

        /// <summary>
        /// The internal array that stores the index into the arrays <see cref="RawValues"/> and <see cref="RawRowIndices"/> of  
        /// the first entry of each column. Its length is equal to <paramref name="NumColumns"/> + 1. 
        /// The last entry is the number of non-zero entries, which must be equal to 
        /// <see cref="RawValues"/>.Length == <see cref="RawRowIndices"/>.Length.
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal int[] RawColOffsets => colOffsets;

        /// <summary>
        /// The internal array that stores the row indices of the non-zero entries in <see cref="RawValues"/>.
        /// Its length is equal to the number of non-zero entries. 
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal int[] RawRowIndices => rowIndices;

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
            => matrixRight.Multiply(vectorLeft, true);
        #endregion

        /// <summary>
        /// See <see cref="IMatrixView.Axpy(IMatrixView, double)"/>.
        /// </summary>
        public IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CscMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    Array.Copy(this.values, resultValues, values.Length);
                    Blas.Daxpy(values.Length, otherCoefficient, otherCSC.values, 0, 1, resultValues, 0, 1);
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
            Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, resultValues, 0, 1);
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
            Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, this.values, 0, 1);
        }

        /// <summary>
        /// See <see cref="IMatrix.Clear"/>.
        /// </summary>
        public void Clear() => Array.Clear(values, 0, values.Length);

        /// <summary>
        /// See <see cref="IMatrixView.Copy(bool)"/>.
        /// </summary>
        IMatrix IMatrixView.Copy(bool copyIndexingData) => Copy(copyIndexingData);

        /// <summary>
        /// Copies the entries of this matrix.
        /// </summary>
        /// <param name="copyIndexingData">
        /// If true, all data of this object will be copied. If false, only the array containing the values of the stored 
        /// matrix entries will be copied. The new matrix will reference the same indexing arrays as this one.
        /// </param>
        public CscMatrix Copy(bool copyIndexingData)
        {
            var valuesCopy = new double[this.values.Length];
            Array.Copy(this.values, valuesCopy, this.values.Length);

            if (!copyIndexingData) return new CscMatrix(NumRows, NumColumns, valuesCopy, this.rowIndices, this.colOffsets);
            else
            {
                var rowIndicesCopy = new int[this.rowIndices.Length];
                Array.Copy(this.rowIndices, rowIndicesCopy, this.rowIndices.Length);
                var colOffsetsCopy = new int[this.colOffsets.Length];
                Array.Copy(this.colOffsets, colOffsetsCopy, this.colOffsets.Length);
                return new CscMatrix(NumRows, NumColumns, valuesCopy, rowIndicesCopy, colOffsetsCopy);
            }
        }

        /// <summary>
        /// See <see cref="IMatrixView.CopyToFullMatrix()"/>
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
        /// See <see cref="IEntrywiseOperableView2D{TMatrixIn, TMatrixOut}.DoEntrywise(TMatrixIn, Func{double, double, double})"/>.
        /// </summary>
        public IMatrix DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
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
        /// See <see cref="IEntrywiseOperable2D{TMatrixIn}.DoEntrywiseIntoThis(TMatrixIn, Func{double, double, double})"/>.
        /// </summary>
        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is CscMatrix casted)
            {
                //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
                if (!HasSameIndexer(casted))
                {
                    throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
                }
                for (int i = 0; i < values.Length; ++i) this.values[i] = binaryOperation(this.values[i], casted.values[i]);
            }
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperableView2D{TMatrixIn, TMatrixOut}.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
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
        /// See <see cref="IEntrywiseOperable2D{TMatrixIn}.DoToAllEntriesIntoThis(Func{double, double})"/>.
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
        /// See <see cref="ISliceable2D.GetColumn(int)"/>.
        /// </summary>
        public Vector GetColumn(int colIndex)
        {
            Preconditions.CheckIndexCol(this, colIndex);
            double[] colVector = new double[NumRows];
            for (int k = colOffsets[colIndex]; k < colOffsets[colIndex + 1]; ++k) colVector[rowIndices[k]] = values[k];
            return Vector.CreateFromArray(colVector, false);
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetRow(int)"/>.
        /// </summary>
        public Vector GetRow(int rowIndex)
        {
            Preconditions.CheckIndexRow(this, rowIndex);
            double[] rowVector = new double[NumColumns];
            for (int j = 0; j < NumColumns; ++j)
            {
                int entryOffset = FindOffsetOf(rowIndex, j);
                if (entryOffset != zeroEntryOffset) rowVector[j] = values[entryOffset];
            }
            return Vector.CreateFromArray(rowVector, false);
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
        /// See <see cref="ISliceable2D.GetSubmatrix(int[], int[])"/>.
        /// </summary>
        public IMatrix GetSubmatrix(int[] rowIndices, int[] colIndices)
            => DenseStrategies.GetSubmatrix(this, rowIndices, colIndices);

        /// <summary>
        /// See <see cref="ISliceable2D.GetSubmatrix(int, int, int, int)"/>.
        /// </summary>
        public IMatrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
            => DenseStrategies.GetSubmatrix(this, rowStartInclusive, rowEndExclusive, colStartInclusive, colEndExclusive);

        /// <summary>
        /// See <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>.
        /// </summary>
        public IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CscMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    if (thisCoefficient == 1.0)
                    {
                        Array.Copy(this.values, resultValues, values.Length);
                        Blas.Daxpy(values.Length, otherCoefficient, otherCSC.values, 0, 1, this.values, 0, 1);
                    }
                    else if (otherCoefficient == 1.0)
                    {
                        Array.Copy(otherCSC.values, resultValues, values.Length);
                        Blas.Daxpy(values.Length, thisCoefficient, this.values, 0, 1, resultValues, 0, 1);
                    }
                    else
                    {
                        Array.Copy(this.values, resultValues, values.Length);
                        BlasExtensions.Daxpby(values.Length, otherCoefficient, otherCSC.values, 0, 1,
                            thisCoefficient, resultValues, 0, 1);
                    }
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
            if (thisCoefficient == 1.0)
            {
                Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, this.values, 0, 1);
            }
            else
            {
                BlasExtensions.Daxpby(values.Length, otherCoefficient, otherMatrix.values, 0, 1, thisCoefficient, this.values, 0, 1);
            }
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyLeft(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            //TODO: To use BLAS for this too, we must accept row major matrices as output.
            if (transposeOther)
            {
                if (transposeThis)
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumColumns);
                    var result = Matrix.CreateZero(other.NumColumns, this.NumRows);
                    CsrMultiplications.MatrixTransTimesCsr(this.NumColumns, values, colOffsets, rowIndices, other, result);
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumRows);
                    var result = Matrix.CreateZero(other.NumColumns, this.NumColumns);
                    CsrMultiplications.MatrixTransTimesCsrTrans(this.NumColumns, values, colOffsets, rowIndices, other, result);
                    return result;
                }
            }
            else
            {
                if (transposeThis)
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumColumns);
                    var result = Matrix.CreateZero(other.NumRows, this.NumRows);
                    CsrMultiplications.MatrixTimesCsr(this.NumColumns, values, colOffsets, rowIndices, other, result);
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumRows);
                    var result = Matrix.CreateZero(other.NumRows, this.NumColumns);
                    CsrMultiplications.MatrixTimesCsrTrans(this.NumColumns, values, colOffsets, rowIndices, other, result);
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
            if (transposeOther)
            {
                if (transposeThis)
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumColumns);
                    var result = Matrix.CreateZero(this.NumColumns, other.NumRows);
                    CsrMultiplications.CsrTimesMatrixTrans(this.NumColumns, values, colOffsets, rowIndices, other, result);
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumColumns);
                    var result = Matrix.CreateZero(this.NumRows, other.NumRows);
                    CsrMultiplications.CsrTransTimesMatrixTrans(this.NumColumns, values, colOffsets, rowIndices, other, result);
                    return result;
                }
            }
            else
            {
                //TODO: perhaps I can use the left multiplications if the other matrix is also transposed
                if (other is Matrix dense) return MultiplyRight(dense, transposeThis);

                if (transposeThis)
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
                    var result = Matrix.CreateZero(this.NumColumns, other.NumColumns);
                    CsrMultiplications.CsrTimesMatrix(this.NumColumns, values, colOffsets, rowIndices, other, result);
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
                    var result = Matrix.CreateZero(this.NumRows, other.NumColumns);
                    CsrMultiplications.CsrTransTimesMatrix(this.NumColumns, values, colOffsets, rowIndices, other, result);
                    return result;
                }
            }
        }

        /// <summary>
        /// Performs the matrix-matrix multiplication: oper(this) * <paramref name="other"/>.
        /// </summary>
        /// <param name="other">
        /// A matrix such that the <see cref="IIndexable2D.NumRows"/> of <paramref name="other"/> are equal to the 
        /// <see cref="IIndexable2D.NumColumns"/> of oper(this).
        /// </param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="otherMatrix"/> has different <see cref="IIndexable2D.NumRows"/> than the 
        /// <see cref="IIndexable2D.NumColumns"/> of oper(this).
        /// </exception>
        public Matrix MultiplyRight(Matrix other, bool transposeThis)
        {
            int numRowsResult;
            if (transposeThis)
            {
                Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
                numRowsResult = this.NumColumns;
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
                numRowsResult = this.NumRows;
            }

            var result = Matrix.CreateZero(numRowsResult, other.NumColumns);
            SparseBlas.Dcscgemm(transposeThis, this.NumRows, other.NumColumns, this.NumColumns, values, colOffsets, rowIndices,
                other.RawData, result.RawData);
            return result;
        }

        /// <summary>
        /// See <see cref="IMatrixView.Multiply(IVectorView, bool)"/>.
        /// </summary>
        public IVector Multiply(IVectorView vector, bool transposeThis = false)
        {
            if (vector is Vector dense) return Multiply(dense, transposeThis);

            if (transposeThis)
            {
                var result = new double[NumColumns];
                Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
                CsrMultiplications.CsrTimesVector(NumColumns, values, colOffsets, rowIndices, vector, result);
                return Vector.CreateFromArray(result, false);
            }
            else
            {
                var result = new double[NumRows];
                Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
                CsrMultiplications.CsrTransTimesVector(NumColumns, values, colOffsets, rowIndices, vector, result);
                return Vector.CreateFromArray(result, false);
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
        public Vector Multiply(Vector vector, bool transposeThis = false)
        {
            //TODO: this performs redundant dimension checks, including checking the transposeThis flag.
            var result = Vector.CreateZero(transposeThis ? NumColumns : NumRows);
            MultiplyIntoResult(vector, result, transposeThis);
            return result;
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyIntoResult(IVectorView, IVector, bool)"/>.
        /// </summary>
        public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
        {
            if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
            {
                MultiplyIntoResult(lhsDense, rhsDense, transposeThis);
            }

            if (transposeThis)
            {
                Preconditions.CheckMultiplicationDimensions(NumRows, lhsVector.Length);
                Preconditions.CheckSystemSolutionDimensions(NumColumns, rhsVector.Length);
                CsrMultiplications.CsrTimesVector(NumColumns, values, colOffsets, rowIndices, lhsVector, rhsVector);
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
                Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
                CsrMultiplications.CsrTransTimesVector(NumColumns, values, colOffsets, rowIndices, lhsVector, rhsVector);
            }
        }

        /// <summary>
        /// Performs the matrix-vector multiplication: <paramref name="rhsVector"/> = oper(this) * <paramref name="vector"/>.
        /// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
        /// The resulting vector will overwrite the entries of <paramref name="rhsVector"/>.
        /// </summary>
        /// <param name="lhsVector">
        /// The vector that will be multiplied by this matrix. It sits on the left hand side of the equation y = oper(A) * x.
        /// Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
        /// == oper(this).<see cref="IIndexable2D.NumColumns"/>.
        /// </param>
        /// <param name="rhsVector">
        /// The vector that will be overwritten by the result of the multiplication. It sits on the right hand side of the 
        /// equation y = oper(A) * x. Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
        /// == oper(this).<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if the <see cref="IIndexable1D.Length"/> of <paramref name="lhsVector"/> or <paramref name="rhsVector"/> 
        /// violate the described contraints.
        /// </exception>
        public void MultiplyIntoResult(Vector lhsVector, Vector rhsVector, bool transposeThis = false)
        {
            if (transposeThis)
            {
                Preconditions.CheckMultiplicationDimensions(NumRows, lhsVector.Length);
                Preconditions.CheckSystemSolutionDimensions(NumColumns, rhsVector.Length);
                
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
                Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
            }
            SparseBlas.Dcscgemv(transposeThis, NumRows, NumColumns, values, colOffsets, rowIndices,
                    lhsVector.RawData, 0, rhsVector.RawData, 0);
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
        IMatrix IMatrixView.Scale(double scalar) => Scale(scalar);

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
            Blas.Dscal(nnz, scalar, resultValues, 0, 1);
            return new CscMatrix(this.NumRows, this.NumColumns, resultValues, this.rowIndices, this.colOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrix.ScaleIntoThis(double)"/>.
        /// </summary>
        public void ScaleIntoThis(double scalar) => Blas.Dscal(values.Length, scalar, values, 0, 1);

        /// <summary>
        /// See <see cref="IMatrix.SetEntryRespectingPattern(int, int, double)"/>.
        /// </summary>
        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            int entryOfsset = FindOffsetOf(rowIdx, colIdx);
            if (entryOfsset == zeroEntryOffset) throw new SparsityPatternModifiedException(
                $"Cannot write to zero entry ({rowIdx}, {colIdx}).");
            else values[entryOfsset] = value;
        }

        /// <summary>
        /// See <see cref="IMatrixView.Transpose"/>.
        /// </summary>
        public IMatrix Transpose() => TransposeToCSR(true);

        /// <summary>
        /// Creates a new <see cref="CscMatrix"/> instance, that is transpose to this: result[i, j] = this[j, i].
        /// </summary>
        public CscMatrix TransposeToCSC()
        {
            // Use C# port of the scipy method.
            // TODO: Perhaps it could be done faster by making extra assumptions. Otherwise use SparseBLAS
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
        /// Otherwise returns <see cref="zeroEntryOffset"/>.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        private int FindOffsetOf(int rowIdx, int colIdx)
        {
            //TODO: if we have a true bool flag to indicate that the row indices of each column are sorted, then use binary search.
            Preconditions.CheckIndices(this, rowIdx, colIdx); //TODO: check indices?
            int colStart = colOffsets[colIdx]; //inclusive
            int colEnd = colOffsets[colIdx + 1]; //exclusive
            for (int k = colStart; k < colEnd; ++k)
            {
                if (rowIndices[k] == rowIdx) return k;
            }
            return zeroEntryOffset;
        }

        private bool HasSameIndexer(CscMatrix other)
        {
            return (this.rowIndices == other.rowIndices) && (this.colOffsets == other.colOffsets);
        }
    }
}
