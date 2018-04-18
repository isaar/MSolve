using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
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
    /// This class is optimized for matrix-vector and matrix-matrix multiplications, where the CSC matrix is on the right 
    /// untransposed or on the left transposed. The other combinations are more efficient using <see cref="CSRMatrix"/>. To build
    /// a <see cref="CSCMatrix"/> conveniently, use <see cref="Builders.DOKColMajor"/>. 
    /// </summary>
    public class CSCMatrix: IMatrix, ISparseMatrix //TODO: Use MKL with descriptors
    {
        public static bool WriteRawArrays { get; set; } = true;

        private readonly double[] values;
        private readonly int[] rowIndices;
        private readonly int[] colOffsets;

        private CSCMatrix(int numRows, int numCols, double[] values, int[] rowIndices, int[] colOffsets)
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

        public double this[int rowIdx, int colIdx]
        {
            get
            {
                int index = FindIndexOf(rowIdx, colIdx);
                if (index == -1) return 0.0;
                else return values[index];
            }
        }

        public static CSCMatrix CreateFromArrays(int numRows, int numCols, double[] values, int[] rowIndices, int[] colOffsets,
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
                if (colOffsets[colOffsets.Length - 1] != values.Length)
                {
                    throw new ArgumentException("The last entry of the CSC column offsets array must be equal to the number of"
                        + " non zero entries, but was " + colOffsets[colOffsets.Length - 1]);
                }
            }
            return new CSCMatrix(numRows, numCols, values, rowIndices, colOffsets);
        }

        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CSCMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    Array.Copy(this.values, resultValues, values.Length);
                    CBlas.Daxpy(values.Length, otherCoefficient, ref otherCSC.values[0], 1, ref resultValues[0], 1);
                    return new CSCMatrix(NumRows, NumColumns, resultValues, this.rowIndices, this.colOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, 1.0, otherMatrix, otherCoefficient);
        }

        public CSCMatrix Axpy(CSCMatrix otherMatrix, double otherCoefficient)
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
            return new CSCMatrix(NumRows, NumColumns, resultValues, this.rowIndices, this.colOffsets);
        }

        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CSCMatrix casted) AxpyIntoThis(casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                 "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        public void AxpyIntoThis(CSCMatrix otherMatrix, double otherCoefficient)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (!HasSameIndexer(otherMatrix))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            CBlas.Daxpy(values.Length, otherCoefficient, ref otherMatrix.values[0], 1, ref this.values[0], 1);
        }

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

        public int CountNonZeros()
        {
            return values.Length;
        }

        public IMatrixView DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is CSCMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherCSC.values[i]);
                    }
                    return new CSCMatrix(NumRows, NumColumns, resultValues, rowIndices, colOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.DoEntrywise(this, other, binaryOperation);
        }

        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is CSCMatrix casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        public void DoEntrywiseIntoThis(CSCMatrix other, Func<double, double, double> binaryOperation)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (!HasSameIndexer(other))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            for (int i = 0; i < values.Length; ++i) this.values[i] = binaryOperation(this.values[i], other.values[i]);
        }

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
                return new CSCMatrix(NumRows, NumColumns, newValues, rowIndicesCopy, colOffsetsCopy);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                return new CSCMatrix(NumRows, NumColumns, newValues, rowIndices, colOffsets).CopyToFullMatrix();
            }
        }

        void IMatrix.DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            DoToAllEntriesIntoThis(unaryOperation);
        }

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

        public SparseFormat GetSparseFormat()
        {
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Row indices", rowIndices);
            format.RawIndexArrays.Add("Column offsets", colOffsets);
            return format;
        }

        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CSCMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    Array.Copy(this.values, resultValues, values.Length);
                    CBlas.Daxpby(values.Length, otherCoefficient, ref otherCSC.values[0], 1, thisCoefficient, ref this.values[0], 1);
                    return new CSCMatrix(NumRows, NumColumns, resultValues, this.rowIndices, this.colOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, thisCoefficient, otherMatrix, otherCoefficient);
        }

        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CSCMatrix casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        public void LinearCombinationIntoThis(double thisCoefficient, CSCMatrix otherMatrix, double otherCoefficient)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (!HasSameIndexer(otherMatrix))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            CBlas.Daxpby(values.Length, otherCoefficient, ref otherMatrix.values[0], 1, thisCoefficient, ref this.values[0], 1);
        }

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

        public VectorMKL MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            if (vector is VectorMKL) return MultiplyRight((VectorMKL)vector, transposeThis);
            else throw new NotImplementedException();
        }

        public VectorMKL MultiplyRight(VectorMKL vector, bool transposeThis = false)
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
                return VectorMKL.CreateFromArray(result, false);
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
                return VectorMKL.CreateFromArray(result, false);
            }
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            int nnz = values.Length;
            for (int i = 0; i < nnz; ++i) aggregator = processEntry(values[i], aggregator);
            aggregator = processZeros(NumRows * NumColumns - nnz, aggregator);
            return finalize(aggregator);
        }

        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            int index = FindIndexOf(rowIdx, colIdx);
            if (index == -1) throw new SparsityPatternModifiedException($"Cannot write to zero entry ({rowIdx}, {colIdx}).");
            else values[index] = value;
        }

        public IMatrixView Transpose()
        {
            return Transpose(true);
        }

        public CSRMatrix Transpose(bool copyInternalArrays)
        {
            if (copyInternalArrays)
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                int[] rowIndicesCopy = new int[rowIndices.Length];
                Array.Copy(rowIndices, rowIndicesCopy, rowIndices.Length);
                int[] colOffsetsCopy = new int[colOffsets.Length];
                Array.Copy(colOffsets, colOffsetsCopy, colOffsets.Length);
                return CSRMatrix.CreateFromArrays(NumColumns, NumRows, valuesCopy, rowIndicesCopy, colOffsetsCopy, false);
            }
            else return CSRMatrix.CreateFromArrays(NumColumns, NumRows, values, rowIndices, colOffsets, false);
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

        private bool HasSameIndexer(CSCMatrix other)
        {
            return (this.rowIndices == other.rowIndices) && (this.colOffsets == other.colOffsets);
        }
    }
}
