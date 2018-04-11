using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;


// TODO: try to make general versions of row major and col major multiplication. The lhs matrix/vector will be supplied by the 
// caller, depending on if it is a matrix, a transposed matrix, a vector, etc. Compare its performance with the verbose code and 
// also try inlining. C preprocessor macros would be actually useful here. Otherwise, move all that boilerplate code to a 
// CSRStrategies static class.
// TODO: In matrix-matrix/vector multiplications: perhaps I should work with a column major array directly instead of an output   
// Matrix and an array instead of an output Vector.
// TODO: perhaps optimizations if (other is Matrix) are needed, to directly index into its raw col major array.
// The access paterns are always the same
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    /// This class is optimized for matrix-vector and matrix-matrix multiplications, where the CSC matrix is on the right 
    /// transposed or on the left untransposed. The other combinations are more efficient using <see cref="CSCMatrix"/>. To build
    /// a <see cref="CSRMatrix"/> conveniently, use <see cref="Builders.DOKRowMajor"/>. 
    /// </summary>
    public class CSRMatrix: IMatrix, ISparseMatrix //TODO: Use MKL with descriptors
    {
        public static bool WriteRawArrays { get; set; } = true;

        private readonly double[] values;
        private readonly int[] colIndices;
        private readonly int[] rowOffsets;

        private CSRMatrix(int numRows, int numCols, double[] values, int[] colIndices, int[] rowOffsets)
        {
            this.values = values;
            this.colIndices = colIndices;
            this.rowOffsets = rowOffsets;
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

        public static CSRMatrix CreateFromArrays(int numRows, int numCols, double[] values, int[] colIndices, int[] rowOffsets, 
            bool checkInput)
        {
            if (checkInput)
            {
                if (rowOffsets.Length != numRows + 1)
                {
                    throw new ArgumentException("The length of the CSR row offsets array must be equal to the number of rows + 1"
                        + ", but was " + rowOffsets.Length);
                }
                if (values.Length != colIndices.Length)
                {
                    throw new ArgumentException("The length of the CSR values and column indices arrays must be equal (and equal"
                        + $" to the number of non zero entries), but were {values.Length} and {colIndices.Length} respectively");
                }
                if (rowOffsets[rowOffsets.Length-1] != values.Length)
                {
                    throw new ArgumentException("The last entry of the CSR row offsets array must be equal to the number of non"
                        + " zero entries, but was " + rowOffsets[rowOffsets.Length - 1]);
                }
            }
            return new CSRMatrix(numRows, numCols, values, colIndices, rowOffsets);
        }

        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CSRMatrix otherCSR) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSR))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    Array.Copy(this.values, resultValues, values.Length);
                    CBlas.Daxpy(values.Length, otherCoefficient, ref otherCSR.values[0], 1, ref resultValues[0], 1);
                    return new CSRMatrix(NumRows, NumColumns, resultValues, this.colIndices, this.rowOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, 1.0, otherMatrix, otherCoefficient);
        }

        public CSRMatrix Axpy(CSRMatrix otherMatrix, double otherCoefficient)
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
            return new CSRMatrix(NumRows, NumColumns, resultValues, this.colIndices, this.rowOffsets);
        }

        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CSRMatrix casted) AxpyIntoThis(casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                 "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        public void AxpyIntoThis(CSRMatrix otherMatrix, double otherCoefficient)
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
            for (int i = 0; i < NumRows; ++i) //Row major order
            {
                int rowStart = rowOffsets[i]; //inclusive
                int rowEnd = rowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    fullMatrix[i, colIndices[k]] = values[k];
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
            if (other is CSRMatrix otherCSR) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSR))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherCSR.values[i]);
                    }
                    return new CSRMatrix(NumRows, NumColumns, resultValues, colIndices, rowOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.DoEntrywise(this, other, binaryOperation);
        }

        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is CSRMatrix casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        public void DoEntrywiseIntoThis(CSRMatrix other, Func<double, double, double> binaryOperation)
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
                int[] colIndicesCopy = new int[colIndices.Length];
                Array.Copy(colIndices, colIndicesCopy, colIndices.Length);
                int[] rowOffsetsCopy = new int[rowOffsets.Length];
                Array.Copy(rowOffsets, rowOffsetsCopy, rowOffsets.Length);
                return new CSRMatrix(NumRows, NumColumns, newValues, colIndicesCopy, rowOffsetsCopy);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                return new CSRMatrix(NumRows, NumColumns, newValues, colIndices, rowOffsets).CopyToFullMatrix();
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
            for (int i = 0; i < NumRows; ++i)
            {
                int rowStart = rowOffsets[i]; //inclusive
                int rowEnd = rowOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    yield return (i, colIndices[k], values[k]);
                }
            }
        }

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
            var comparer = new ValueComparer(1e-13);
            for (int i = 0; i < NumRows; ++i)
            {
                int rowStart = rowOffsets[i]; //inclusive
                int rowEnd = rowOffsets[i + 1]; //exclusive
                int previousCol = 0;
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    int col = colIndices[k];
                    for (int j = previousCol; j < col; ++j) // zero entries between the stored ones
                    {
                        if (!comparer.AreEqual(0.0, other[i, j])) return false;
                    }
                    if (!comparer.AreEqual(values[k], other[i, col])) return false; // Non zero entry
                    previousCol = col + 1;
                }
            }
            return true; //At this point all entries have been checked and are equal
        }

        public SparseFormat GetSparseFormat()
        {
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Column indices", colIndices);
            format.RawIndexArrays.Add("Row offsets", rowOffsets);
            return format;
        }

        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CSRMatrix otherCSR) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherCSR))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    Array.Copy(this.values, resultValues, values.Length);
                    CBlas.Daxpby(values.Length, otherCoefficient, ref otherCSR.values[0], 1, thisCoefficient, ref this.values[0], 1);
                    return new CSRMatrix(NumRows, NumColumns, resultValues, this.colIndices, this.rowOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, thisCoefficient, otherMatrix, otherCoefficient);
        }

        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is CSRMatrix casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        public void LinearCombinationIntoThis(double thisCoefficient, CSRMatrix otherMatrix, double otherCoefficient)
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
            // TODO: Throwing exceptions when csr is on the right seems attractive.
            if (transposeThis) // Compute one output column at a time. TODO: perhaps 1 row at a time.
            {
                if (transposeOther)
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumColumns);
                    var result = Matrix.CreateZero(other.NumColumns, this.NumRows);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
                    {
                        int rowStart = rowOffsets[c]; //inclusive
                        int rowEnd = rowOffsets[c + 1]; //exclusive
                        for (int i = 0; i < other.NumColumns; ++i)
                        {
                            double dot = 0.0;
                            for (int k = rowStart; k < rowEnd; ++k) // other.col[i] * transpose(csr).col[j] = other.col[i] * csr.row[j]
                            {
                                dot += values[k] * other[colIndices[k], i]; 
                            }
                            result[i, c] = dot;
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumColumns);
                    var result = Matrix.CreateZero(other.NumRows, this.NumRows);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time.
                    {
                        int csrRowStart = rowOffsets[c]; //inclusive
                        int csrRowEnd = rowOffsets[c + 1]; //exclusive
                        for (int i = 0; i < other.NumRows; ++i)
                        {
                            double dot = 0.0;
                            for (int k = csrRowStart; k < csrRowEnd; ++k) // other.row[i] * transpose(csr).col[j] = other.row[i] * csr.row[j]
                            {
                                dot += values[k] * other[i, colIndices[k]]; 
                            }
                            result[i, c] = dot;
                        }
                    }
                    return result;
                }
            }
            else // Compute one output row at a time. TODO: not efficent for the col major output matrix.
            {
                if (transposeOther)
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumRows, this.NumRows);
                    var result = Matrix.CreateZero(other.NumColumns, this.NumColumns);
                    for (int r = 0; r < result.NumRows; ++r) // Compute one output row at a time.
                    {
                        // x * A = linear combination of rows of A with the entries of x as coefficients,
                        // where x is row r of transpose(other matrix).
                        for (int i = 0; i < this.NumRows; ++i)
                        {
                            double scalar = other[i, r];
                            int csrRowStart = rowOffsets[i]; //inclusive
                            int csrRowEnd = rowOffsets[i + 1]; //exclusive
                            for (int k = csrRowStart; k < csrRowEnd; ++k) // sum(other[i,r] * csr.row[i]))
                            {
                                result[r, colIndices[k]] += scalar * values[k];
                            }
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumRows);
                    var result = Matrix.CreateZero(other.NumRows, this.NumColumns);
                    for (int r = 0; r < result.NumRows; ++r) // Compute one output row at a time.
                    {
                        // x * A = linear combination of rows of A with the entries of x as coefficients,
                        // where x is row r of the other matrix.
                        for (int i = 0; i < this.NumRows; ++i)
                        {
                            double scalar = other[r, i];
                            int csrRowStart = rowOffsets[i]; //inclusive
                            int csrRowEnd = rowOffsets[i + 1]; //exclusive
                            for (int k = csrRowStart; k < csrRowEnd; ++k) // sum(other[r,i] * csr.row[i]))
                            {
                                result[r, colIndices[k]] += scalar * values[k];
                            }
                        }
                    }
                    return result;
                }
            }
        }

        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            if (transposeThis)
            {
                if (transposeOther)
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumColumns);
                    var result = Matrix.CreateZero(this.NumColumns, other.NumRows);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
                    {
                        // A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients, 
                        // where x is column c of the other matrix
                        for (int i = 0; i < this.NumRows; ++i)
                        {
                            double scalar = other[c, i]; 
                            int csrRowStart = rowOffsets[i]; //inclusive
                            int csrRowEnd = rowOffsets[i + 1]; //exclusive
                            for (int k = csrRowStart; k < csrRowEnd; ++k) // sum(other[c,i] * transpose(csr.col[c])) = sum(other[c,i] * csr.row[c])
                            {
                                result[colIndices[k], c] += scalar * values[k];
                            }
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
                    var result = Matrix.CreateZero(this.NumColumns, other.NumColumns);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
                    {
                        // A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients, 
                        // where x is column c of the other matrix
                        for (int i = 0; i < this.NumRows; ++i)
                        {
                            double scalar = other[i, c]; 
                            int csrRowStart = rowOffsets[i]; //inclusive
                            int csrRowEnd = rowOffsets[i + 1]; //exclusive
                            for (int k = csrRowStart; k < csrRowEnd; ++k) // sum(other[i,c] * transpose(csr.col[c])) = sum(other[i,c] * csr.row[c])
                            {
                                result[colIndices[k], c] += scalar * values[k];
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
                    Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumColumns);
                    var result = Matrix.CreateZero(this.NumRows, other.NumRows);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
                    {
                        for (int i = 0; i < this.NumRows; ++i)
                        {
                            double dot = 0.0;
                            int csrRowStart = rowOffsets[i]; //inclusive
                            int csrRowEnd = rowOffsets[i + 1]; //exclusive
                            for (int k = csrRowStart; k < csrRowEnd; ++k) // csr.row[i] * other.row[c]
                            {
                                dot += values[k] * other[c, colIndices[k]];
                            }
                            result[i, c] = dot;
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
                    var result = Matrix.CreateZero(this.NumRows, other.NumColumns);
                    for (int c = 0; c < result.NumColumns; ++c) // Compute one output column at a time
                    {
                        for (int i = 0; i < this.NumRows; ++i)
                        {
                            double dot = 0.0;
                            int csrRowStart = rowOffsets[i]; //inclusive
                            int csrRowEnd = rowOffsets[i + 1]; //exclusive
                            for (int k = csrRowStart; k < csrRowEnd; ++k) // csr.row[i] * other.col[c]
                            {
                                dot += values[k] * other[colIndices[k], c];
                            }
                            result[i, c] = dot;
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

        /// <summary>
        /// 
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="transposeThis">If this method is called multiple times with transposeThis = true, consider using a CSC instead</param>
        /// <returns></returns>
        public VectorMKL MultiplyRight(VectorMKL vector, bool transposeThis = false)
        {
            if (transposeThis)
            {
                Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
                // A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
                double[] result = new double[NumColumns];
                for (int i = 0; i < NumRows; ++i)
                {
                    double scalar = vector[i];
                    int rowStart = rowOffsets[i]; //inclusive
                    int rowEnd = rowOffsets[i + 1]; //exclusive
                    for (int k = rowStart; k < rowEnd; ++k)
                    {
                        result[colIndices[k]] += scalar * values[k];
                    }
                }
                return VectorMKL.CreateFromArray(result, false);
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
                double[] result = new double[NumRows];
                for (int i = 0; i < NumRows; ++i)
                {
                    double dot = 0.0;
                    int rowStart = rowOffsets[i]; //inclusive
                    int rowEnd = rowOffsets[i + 1]; //exclusive
                    for (int k = rowStart; k < rowEnd; ++k)
                    {
                        dot += values[k] * vector[colIndices[k]];
                    }
                    result[i] = dot;
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

        public CSCMatrix Transpose(bool copyInternalArrays )
        {
            if (copyInternalArrays)
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                int[] colIndicesCopy = new int[colIndices.Length];
                Array.Copy(colIndices, colIndicesCopy, colIndices.Length);
                int[] rowOffsetsCopy = new int[rowOffsets.Length];
                Array.Copy(rowOffsets, rowOffsetsCopy, rowOffsets.Length);
                return CSCMatrix.CreateFromArrays(NumColumns, NumRows, valuesCopy, colIndicesCopy, rowOffsetsCopy, false);
            }
            else return CSCMatrix.CreateFromArrays(NumColumns, NumRows, values, colIndices, rowOffsets, false);
        }

        /// <summary>
        /// Return the index into values and colIndices arrays, if the (rowIdx, colIdx) entry is within the pattern. 
        /// Otherwise returns -1.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        /// <returns></returns>
        private int FindIndexOf(int rowIdx, int colIdx)
        {
            Preconditions.CheckIndices(this, rowIdx, colIdx); //TODO: check indices?
            int rowStart = rowOffsets[rowIdx]; //inclusive
            int rowEnd = rowOffsets[rowIdx + 1]; //exclusive
            for (int k = rowStart; k < rowEnd; ++k)
            {
                if (colIndices[k] == colIdx) return k;
            }
            return -1;
        }

        private bool HasSameIndexer(CSRMatrix other)
        {
            return (this.colIndices == other.colIndices) && (this.rowOffsets == other.rowOffsets);
        }
    }
}
