using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    //TODO: Use MKL with descriptors
    public class CSCMatrix: IEntrywiseOperable, IIndexable2D, IWriteable
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

        public IEntrywiseOperable DoEntrywise(IEntrywiseOperable other, Func<double, double, double> binaryOperation)
        {
            if (other is CSCMatrix otherCSC) // In case both matrices have the exact same index arrays
            {
                if ((this.rowIndices == otherCSC.rowIndices) && (this.colOffsets == otherCSC.colOffsets))
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

        public void DoEntrywiseIntoThis(CSCMatrix other, Func<double, double, double> binaryOperation)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if ((this.rowIndices != other.rowIndices) || (this.colOffsets != other.colOffsets))
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

        public void WriteToConsole()
        {
            if (WriteRawArrays)
            {
                Console.WriteLine("Values:");
                for (int i = 0; i < values.Length; ++i)
                {
                    Console.Write(values[i]);
                    Console.Write(" ");
                }
                Console.WriteLine("Row indices:");
                for (int i = 0; i < rowIndices.Length; ++i)
                {
                    Console.Write(rowIndices[i]);
                    Console.Write(" ");
                }
                Console.WriteLine("Column offsets:");
                for (int i = 0; i < colOffsets.Length; ++i)
                {
                    Console.Write(colOffsets[i]);
                    Console.Write(" ");
                }
            }
            else
            {
                var formatter = new SparseMatrixFormatting();
                for (int j = 0; j < NumColumns; ++j)
                {
                    int colStart = colOffsets[j]; //inclusive
                    int colEnd = colOffsets[j + 1]; //exclusive
                    for (int k = colStart; k < colEnd; ++k)
                    {
                        Console.WriteLine(formatter.FormatNonZeroEntry(rowIndices[k], j, values[k]));
                    }
                }
            }
        }

        /// <summary>
        /// Write the entries of the matrix to a specified file. If the file doesn't exist a new one will be created.
        /// </summary>
        /// <param name="path">The path of the file and its extension.</param>
        /// <param name="append">If the file already exists: Pass <see cref="append"/> = true to write after the current end of 
        ///     the file. Pass<see cref="append"/> = false to overwrite the file.</param>
        public void WriteToFile(string path, bool append = false)
        {
            //TODO: incorporate this and WriteToConsole into a common function, where the user passes the stream and an object to 
            //deal with formating. Also add support for relative paths. Actually these methods belong in the "Logging" project, 
            // but since they are extremely useful they are implemented here for now.
            using (var writer = new StreamWriter(path, append))
            {
                if (WriteRawArrays)
                {
                    writer.WriteLine("Values:");
                    for (int i = 0; i < values.Length; ++i)
                    {
                        writer.Write(values[i]);
                        writer.Write(" ");
                    }
                    writer.WriteLine("Row indices:");
                    for (int i = 0; i < rowIndices.Length; ++i)
                    {
                        writer.Write(rowIndices[i]);
                        writer.Write(" ");
                    }
                    writer.WriteLine("Column offsets:");
                    for (int i = 0; i < colOffsets.Length; ++i)
                    {
                        writer.Write(colOffsets[i]);
                        writer.Write(" ");
                    }
                }
                else
                {
                    var formatter = new SparseMatrixFormatting();
                    for (int j = 0; j < NumColumns; ++j)
                    {
                        int colStart = colOffsets[j]; //inclusive
                        int colEnd = colOffsets[j + 1]; //exclusive
                        for (int k = colStart; k < colEnd; ++k)
                        {
                            writer.WriteLine(formatter.FormatNonZeroEntry(rowIndices[k], j, values[k]));
                        }
                    }
                }

#if DEBUG
                writer.Flush(); // If the user inspects the file while debugging, make sure the contentss are written.
#endif
            }
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
    }
}
