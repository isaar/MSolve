using System;
using System.Runtime.CompilerServices;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using static ISAAR.MSolve.LinearAlgebra.LibrarySettings;

//TODO: Should this be removed? I need it to provide analyzers with a matrix.
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Symmetric sparse matrix in Compressed Sparse Columns format, with only the non-zero entries of the upper triangle being 
    /// explicitly stored. This matrix format is better used for factorizations using SuiteSparse or CSparse libraries.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SymmetricCscMatrix: IMatrix
    {
        /// <summary>
        /// values.Length = number of non zeros in upper triangle.
        /// </summary>
        private readonly double[] values;

        /// <summary>
        /// rowIndices.Length = number of non zeros in upper triangle.
        /// </summary>
        private readonly int[] rowIndices;

        /// <summary>
        /// colOffsets.Length = number of rows/columns +1. 
        /// colOffsets[colOffsets.Length-1] = number of non zeros in upper triangle.
        /// </summary>
        private readonly int[] colOffsets;

        private SymmetricCscMatrix(int order, int numNonZerosUpper, double[] values, int[] rowIndices, int[] colOffsets)
        {
            this.NumColumns = order;
            this.NumRows = order;
            this.NumNonZerosUpper = numNonZerosUpper;
            this.values = values;
            this.rowIndices = rowIndices;
            this.colOffsets = colOffsets;
        }

        /// <summary>
        /// The internal array that stores the non-zero entries of the upper triangle. The non-zero entries of each 
        /// column are consecutive. Its length is equal to the number of the upper triangle's non-zero entries. 
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal double[] RawValues => values;

        /// <summary>
        /// The internal array that stores the index into the arrays <see cref="RawValues"/> and <see cref="RawRowIndices"/> of  
        /// the first entry of each column. Its length is equal to <paramref name="NumColumns"/> + 1. 
        /// The last entry is the number of the upper triangle's non-zero entries, which must be equal to 
        /// <see cref="RawValues"/>.Length == <see cref="RawRowIndices"/>.Length.
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal int[] RawColOffsets => colOffsets;

        /// <summary>
        /// The internal array that stores the row indices of the non-zero entries in <see cref="RawValues"/>.
        /// Its length is equal to the number of the upper triangle's non-zero entries. 
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal int[] RawRowIndices => rowIndices;

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of the upper triangle's non zero entries. These are the only ones being explicitly stored.
        /// </summary>
        public int NumNonZerosUpper { get; }

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
                int offset = FindOffsetOf(rowIdx, colIdx);
                if (offset >= 0) return values[offset];
                else return 0.0;
            }
        }

        /// <summary>
        /// Initializes a new <see cref="SymmetricCscMatrix"/> with the specified dimensions and the provided arrays 
        /// (<paramref name="values"/>, <paramref name="rowIndices"/> and <paramref name="colOffsets"/>) as its internal data.
        /// </summary>
        /// <param name="order">The number of rows /columns of the new matrix.</param>
        /// <param name="values">
        /// Array that contains the non-zero entries of the upper triangle. It must have the same length as 
        /// <paramref name="rowIndices"/>. The non-zero entries of each column must appear consecutively in 
        /// <paramref name="values"/>. They can also be sorted in increasing order of their row indices, which speeds up 
        /// subsequent operations.
        /// </param>
        /// <param name="rowIndices">
        /// Array that contains the row indices of the upper triangle's non-zero entries. It must have the same length as 
        /// <paramref name="values"/>. There is an 1 to 1 matching between these two arrays: <paramref name="rowIndices"/>[i] 
        /// is the row index of the entry <paramref name="values"/>[i]. Also: 
        /// 0 &lt;= <paramref name="rowIndices"/>[i] &lt; <paramref name="numRows"/>.
        /// </param>
        /// <param name="colOffsets">
        /// Array that contains the index of the first entry of each column into the arrays <paramref name="values"/> and 
        /// <paramref name="rowIndices"/>. Its length is <paramref name="numRows"/> + 1. The last entry is the number of 
        /// non-zero entries, which must be equal to the length of <paramref name="values"/>and <paramref name="rowIndices"/>.
        /// </param>
        /// <param name="checkInput">
        /// If true, the provided arrays will be checked to make sure they are valid symmetric CSC arrays, which is safer. 
        /// If false, no such check will take place, which is faster.
        /// </param>
        public static SymmetricCscMatrix CreateFromArrays(int order, double[] values, int[] rowIndices, int[] colOffsets,
            bool checkInput)
        {
            int nnz = colOffsets[colOffsets.Length - 1];
            if (checkInput)
            {
                if (colOffsets.Length != order + 1)
                {
                    throw new ArgumentException("The length of the symmetric CSC column offsets array must be equal to the order"
                        + " of the matrix + 1, but was " + colOffsets.Length);
                }
                if ((nnz != values.Length) || (nnz != rowIndices.Length))
                {
                    throw new ArgumentException("Mismatch in dimensions of the symmetric CSC arrays. Check that"
                        + " colOffsets.Length = number of rows/columns + 1 and that colOffsets[colOffsets.Length-1]"
                        + " = values.Length = rowIndices.Length = number of non zeros in upper triangle");
                }
            }
            return new SymmetricCscMatrix(order, nnz, values, rowIndices, colOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Axpy(IMatrixView, double)"/>.
        /// </summary>
        public IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if ((otherMatrix is SymmetricCscMatrix otherCSC) && HaveSameIndexArrays(otherCSC))
            {
                // Unneeded if the indexers are identical, but worth it
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix); 

                //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
                double[] resultValues = new double[values.Length];
                Array.Copy(this.values, resultValues, values.Length);

                Blas.Daxpy(values.Length, otherCoefficient, otherCSC.values, 0, 1, resultValues, 0, 1);

                // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                return new SymmetricCscMatrix(NumRows, NumColumns, resultValues, this.rowIndices, this.colOffsets);
            }
            else return DoEntrywise(otherMatrix, (thisEntry, otherEntry) => thisEntry + otherCoefficient * otherEntry);
        }

        /// <summary>
        /// See <see cref="IMatrix.AxpyIntoThis(IMatrixView, double)"/>.
        /// </summary>
        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            if ((otherMatrix is SymmetricCscMatrix otherCSC) && HaveSameIndexArrays(otherCSC))
            {
                Blas.Daxpy(values.Length, otherCoefficient, otherCSC.values, 0, 1, this.values, 0, 1);
            }
            else throw new SparsityPatternModifiedException("Only allowed if the index arrays are the same");
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
        public SymmetricCscMatrix Copy(bool copyIndexingArrays)
        {
            var valuesCopy = new double[values.Length];
            Array.Copy(values, valuesCopy, values.Length);
            if (!copyIndexingArrays)
            {
                return new SymmetricCscMatrix(NumColumns, NumNonZerosUpper, valuesCopy, rowIndices, colOffsets);
            }
            else
            {
                var rowIndicesCopy = new int[rowIndices.Length];
                Array.Copy(rowIndices, rowIndicesCopy, rowIndices.Length);
                var colOffsetsCopy = new int[colOffsets.Length];
                Array.Copy(colOffsets, colOffsetsCopy, colOffsets.Length);
                return new SymmetricCscMatrix(NumColumns, NumNonZerosUpper, valuesCopy, rowIndicesCopy, colOffsetsCopy);
            }
        }

        /// <summary>
        /// See <see cref="IMatrixView.CopyToFullMatrix()"/>
        /// </summary>
        public Matrix CopyToFullMatrix()
        {
            Matrix fullMatrix = Matrix.CreateZero(this.NumRows, this.NumColumns);
            for (int j = 0; j < this.NumColumns; ++j) //Column major order
            {
                int colCurrent = colOffsets[j];
                int colNext = colOffsets[j + 1]; // TODO: 1 of the two accesses can be removed
                for (int i = colCurrent; i < colNext; ++i)
                {
                    fullMatrix[rowIndices[i], colCurrent] = values[i];
                }
            }
            return fullMatrix;
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperableView2D{TMatrixIn, TMatrixOut}.DoEntrywise(TMatrixIn, Func{double, double, double})"/>.
        /// </summary>
        public IMatrix DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            if (other is SymmetricCscMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if (HaveSameIndexArrays(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    var resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherCSC.values[i]);
                    }
                    return new SymmetricCscMatrix(NumColumns, NumNonZerosUpper, resultValues, rowIndices, colOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            Matrix result = Matrix.CreateZero(this.NumRows, this.NumColumns);
            for (int j = 0; j < this.NumColumns; ++j)
            {
                for (int i = 0; i < this.NumRows; ++i)
                {
                    result[i, j] = binaryOperation(this[i, j], other[i, j]);
                }
            }
            return result;
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperable2D{TMatrixIn}.DoEntrywiseIntoThis(TMatrixIn, Func{double, double, double})"/>.
        /// </summary>
        public void DoEntrywiseIntoThis(IMatrixView matrix, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, matrix); // Unneeded if they have the same indexers, but worth it
            if ((matrix is SymmetricCscMatrix otherCSC) && HaveSameIndexArrays(otherCSC))
            {
                for (int i = 0; i < values.Length; ++i) this.values[i] = binaryOperation(this.values[i], otherCSC.values[i]);
            }
            else throw new SparsityPatternModifiedException("Only allowed if the index arrays are the same");
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperableView2D{TMatrixIn, TMatrixOut}.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Only apply the operation on non zero entries
                double[] newValues = new double[values.Length];
                for (int i = 0; i < values.Length; ++i) newValues[i] = unaryOperation(values[i]);

                //TODO: Perhaps I should also copy the indexers
                return new SymmetricCscMatrix(NumColumns, NumNonZerosUpper, newValues, rowIndices, colOffsets);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                Matrix full = CopyToFullMatrix();
                full.DoToAllEntriesIntoThis(unaryOperation);
                return full;
            }
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperable2D{TMatrixIn}.DoToAllEntriesIntoThis(Func{double, double})"/>.
        /// </summary>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                for (int i = 0; i < values.Length; ++i) values[i] = unaryOperation(values[i]);
            }
            else
            {
                throw new SparsityPatternModifiedException("This operation will change the sparsity pattern");
            }
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
            => DenseStrategies.AreEqual(this, other, tolerance);

        /// <summary>
        /// See <see cref="ISliceable2D.GetColumn(int)"/>.
        /// </summary>
        public Vector GetColumn(int colIndex) => DenseStrategies.GetColumn(this, colIndex);

        /// <summary>
        /// See <see cref="ISliceable2D.GetRow(int)"/>.
        /// </summary>
        public Vector GetRow(int rowIndex) => GetColumn(rowIndex);

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
            if ((otherMatrix is SymmetricCscMatrix otherCSC) && HaveSameIndexArrays(otherCSC))
            {
                // Unneeded if the indexers are identical, but worth it
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);

                //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
                double[] resultValues = new double[values.Length];

                if (thisCoefficient == 1.0)
                {
                    Array.Copy(this.values, resultValues, values.Length);
                    Blas.Daxpy(values.Length, otherCoefficient, otherCSC.values, 0, 1, resultValues, 0, 1);
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

                // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                return new SymmetricCscMatrix(NumRows, NumColumns, resultValues, this.rowIndices, this.colOffsets);
            }
            else return DoEntrywise(otherMatrix, 
                (thisEntry, otherEntry) => thisCoefficient * thisEntry + otherCoefficient * otherEntry);
        }

        /// <summary>
        /// See <see cref="IMatrix.LinearCombinationIntoThis(double, IMatrixView, double)"/>.
        /// </summary>
        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            if ((otherMatrix is SymmetricCscMatrix otherCSC) && HaveSameIndexArrays(otherCSC))
            {
                Blas.Daxpy(values.Length, otherCoefficient, otherCSC.values, 0, 1, this.values, 0, 1);
            }
            else throw new SparsityPatternModifiedException("Only allowed if the index arrays are the same");
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyLeft(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
            => DenseStrategies.Multiply(other, this, transposeOther, transposeThis);

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
            => DenseStrategies.Multiply(this, other, transposeThis, transposeOther);

        /// <summary>
        /// See <see cref="IMatrixView.Multiply(IVectorView, bool)"/>.
        /// </summary>
        public IVector Multiply(IVectorView vector, bool transposeThis = false)
            => DenseStrategies.Multiply(this, vector, transposeThis);

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyIntoResult(IVectorView, IVector, bool)"/>.
        /// </summary>
        public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis)
            => DenseStrategies.MultiplyIntoResult(this, lhsVector, rhsVector, transposeThis);

        /// <summary>
        /// Matrix-vector multiplication, with the vector on the right: matrix * vector or transpose(matrix) * vector.
        /// </summary>
        /// <param name="vector">A vector with length equal to <see cref="NumColumns"/>.</param>
        /// <param name="transposeThis">Set to true to transpose this (the left matrix). Unless the transpose matrix is used in 
        ///     more than one multiplications, setting this flag to true is usually preferable to creating the transpose.</param>
        /// <returns></returns>
        public Vector MultiplyRight(Vector vector, bool transposeThis = false)
            => DenseStrategies.Multiply(this, vector, transposeThis);

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            int numNonZeros = 0;
            for (int j = 0; j < NumColumns; ++j)
            {
                int colStart = colOffsets[j]; //inclusive
                int colEnd = colOffsets[j + 1]; //exclusive
                for (int k = colStart; k < colEnd; ++k)
                {
                    
                    if (rowIndices[k] == j) 
                    {
                        aggregator = processEntry(values[k], aggregator);
                        ++numNonZeros;
                    }
                    else // Do the above twice for entries not on the diagonal
                    {
                        aggregator = processEntry(values[k], aggregator);
                        aggregator = processEntry(values[k], aggregator);
                        numNonZeros += 2;
                    }
                }
            }
            aggregator = processZeros(NumRows * NumColumns - numNonZeros, aggregator);
            return finalize(aggregator);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Scale(double)"/>.
        /// </summary>
        public IMatrix Scale(double scalar)
        {
            // Only apply the operation on non zero entries
            var resultValues = new double[values.Length];
            Blas.Dscal(NumNonZerosUpper, scalar, resultValues, 0, 1);

            //TODO: Perhaps I should also copy the indexers
            return new SymmetricCscMatrix(NumColumns, NumNonZerosUpper, resultValues, rowIndices, colOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrix.ScaleIntoThis(double)"/>.
        /// </summary>
        public void ScaleIntoThis(double scalar) => Blas.Dscal(NumNonZerosUpper, scalar, this.values, 0, 1);

        /// <summary>
        /// See <see cref="IMatrix.SetEntryRespectingPattern(int, int, double)"/>.
        /// </summary>
        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            int offset = FindOffsetOf(rowIdx, colIdx);
            if (offset >= 0) values[offset] = value;
            else throw new SparsityPatternModifiedException($"Cannot write to zero entry ({rowIdx}, {colIdx}).");
        }

        /// <summary>
        /// See <see cref="IMatrixView.Transpose"/>.
        /// </summary>
        public IMatrix Transpose() => Copy(true);

        /// <summary>
        /// 2 or more sparse matrices may have the same (2 references to 1 object) index arrays. 
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private bool HaveSameIndexArrays(SymmetricCscMatrix other)
            => (this.rowIndices == other.rowIndices) && (this.colOffsets == other.colOffsets);

        /// <summary>
        /// For an entry (i,j), returns the offset into values and rowIndices arrays or -1 of the entry corresponds to a 
        /// structural zero.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        private int FindOffsetOf(int rowIdx, int colIdx)
        {
            if (rowIdx > colIdx)
            {
                int swap = rowIdx;
                rowIdx = colIdx;
                colIdx = swap;
            }
            int colStart = colOffsets[colIdx];
            int colEnd = colOffsets[colIdx + 1];
            for (int k = colStart; k < colEnd; ++k) //Only scan the nnz entries for the given column
            {
                if (rowIndices[k] == rowIdx) return k;
            }
            return -1;
        }
    }
}
