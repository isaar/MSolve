using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    //TODO: Use MKL with descriptors
    public class CSRMatrix: IEntrywiseOperable, IIndexable2D, ISparseMatrix, ITransposable
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

        public IEntrywiseOperable DoEntrywise(IEntrywiseOperable other, Func<double, double, double> binaryOperation)
        {
            if (other is CSRMatrix otherCSR) // In case both matrices have the exact same index arrays
            {
                if ((this.colIndices == otherCSR.colIndices) && (this.rowOffsets == otherCSR.rowOffsets))
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

        public void DoEntrywiseIntoThis(CSRMatrix other, Func<double, double, double> binaryOperation)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if ((this.colIndices != other.colIndices) || (this.rowOffsets != other.rowOffsets))
            {
                throw new SparsityPatternModifiedException("Only allowed if the index arrays are the same");
            }
            for (int i = 0; i < values.Length; ++i) this.values[i] = binaryOperation(this.values[i], other.values[i]);
        }

        IEntrywiseOperable IEntrywiseOperable.DoToAllEntries(Func<double, double> unaryOperation)
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

        public ITransposable Transpose()
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

        
    }
}
