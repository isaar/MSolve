using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    //TODO: Use MKL with descriptors
    public class CSCMatrix: IIndexable2D
    {
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
            this.NumNonZeros = values.Length;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// Only structural non zeros
        /// </summary>
        /// <returns></returns>
        public int NumNonZeros { get; }

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
