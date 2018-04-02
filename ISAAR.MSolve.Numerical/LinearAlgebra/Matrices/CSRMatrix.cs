using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.LinearAlgebra.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    //TODO: Use MKL with descriptors
    public class CSRMatrix: IIndexable2D, IWriteable, ISparseMatrix
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

        public int CountNonZeros()
        {
            return values.Length;
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
                Console.WriteLine("Column indices:");
                for (int i = 0; i < colIndices.Length; ++i)
                {
                    Console.Write(colIndices[i]);
                    Console.Write(" ");
                }
                Console.WriteLine("Row offsets:");
                for (int i = 0; i < rowOffsets.Length; ++i)
                {
                    Console.Write(rowOffsets[i]);
                    Console.Write(" ");
                }
            }
            else
            {
                var formatter = new SparseMatrixFormatting();
                for (int i = 0; i < NumRows; ++i)
                {
                    int rowStart = rowOffsets[i]; //inclusive
                    int rowEnd = rowOffsets[i + 1]; //exclusive
                    for (int k = rowStart; k < rowEnd; ++k)
                    {
                       Console.WriteLine(formatter.FormatNonZeroEntry(i, colIndices[k], values[k]));
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
                    writer.WriteLine("Column indices:");
                    for (int i = 0; i < colIndices.Length; ++i)
                    {
                        writer.Write(colIndices[i]);
                        writer.Write(" ");
                    }
                    writer.WriteLine("Row offsets:");
                    for (int i = 0; i < rowOffsets.Length; ++i)
                    {
                        writer.Write(rowOffsets[i]);
                        writer.Write(" ");
                    }
                }
                else
                {
                    var formatter = new SparseMatrixFormatting();
                    for (int i = 0; i < NumRows; ++i)
                    {
                        int rowStart = rowOffsets[i]; //inclusive
                        int rowEnd = rowOffsets[i + 1]; //exclusive
                        for (int k = rowStart; k < rowEnd; ++k)
                        {
                            writer.WriteLine(formatter.FormatNonZeroEntry(i, colIndices[k], values[k]));
                        }
                    }
                }

#if DEBUG
                writer.Flush(); // If the user inspects the file while debugging, make sure the contentss are written.
#endif
            }
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
