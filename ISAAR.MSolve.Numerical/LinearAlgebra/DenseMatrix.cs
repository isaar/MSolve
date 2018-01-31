using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class DenseMatrix: IMatrix
    {
        protected readonly double[,] data;

        public double this[int row, int col]
        {
            get { return data[row, col]; }
            set { data[row, col] = value; }
        }

        public int Rows { get; }
        public int Columns { get; }
        public int NumNonZeros { get { return Rows * Columns; } }

        #region construction
        protected DenseMatrix(double[,] data)
        {
            this.data = data;
            Rows = data.GetLength(0);
            Columns = data.GetLength(1);
        }

        public static DenseMatrix CreateZero(int rows, int columns)
        {
            return new DenseMatrix(new double[rows, columns]);
        }

        public static DenseMatrix CreateWithValue(int rows, int columns, double value)
        {
            double[,] data = new double[rows, columns];
            for (int r = 0; r < rows; ++r)
            {
                for (int c = 0; c < columns; ++c)
                {
                    data[r, c] = value;
                }
            }
            return new DenseMatrix(data);
        }

        public static DenseMatrix CreateFromArray(double[,] data, bool copyArray = true)
        {
            if (copyArray)
            {
                double[,] clone = new double[data.GetLength(0), data.GetLength(1)];
                for (int r = 0; r < data.GetLength(0); ++r)
                {
                    for (int c = 0; c < data.GetLength(1); ++c)
                    {
                        clone[r, c] = data[r, c];
                    }
                }
                return new DenseMatrix(clone);
            }
            else return new DenseMatrix(data);
        }

        public static DenseMatrix CreateFromMatrix(DenseMatrix matrix)
        {
            return DenseMatrix.CreateFromArray(matrix.data, true);
        }
        #endregion

        public double[,] CopyToArray()
        {
            var clone = new double[data.GetLength(0), data.GetLength(1)];
            Array.Copy(data, clone, data.Length); // 2D Array.Copy works fine, if the destination and source have the same size.
            return clone;
        }

        public IMatrix DoPointwise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckMatrixDimensionsSame(this, other);
            var result = new double[this.Rows, this.Columns];
            for (int r = 0; r < this.Rows; ++r)
            {
                for (int c = 0; c < this.Columns; ++c)
                {
                    result[r, c] = binaryOperation(this.data[r, c], other[r, c]);
                }
            }
            return new DenseMatrix(result);
        }

        public void DoPointwiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckMatrixDimensionsSame(this, other);
            for (int r = 0; r < this.Rows; ++r)
            {
                for (int c = 0; c < this.Columns; ++c)
                {
                    data[r, c] = binaryOperation(this.data[r, c], other[r, c]);
                }
            }
        }

        public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
        {
            var result = new double[Rows, Columns];
            for (int r = 0; r < Rows; ++r)
            {
                for (int c = 0; c < Columns; ++c)
                {
                    result[r, c] = unaryOperation(data[r, c]);
                }
            }
            return new DenseMatrix(result);
        }

        // Ok for a DenseMatrix, but for sparse formats some operation (e.g scale) maintain the sparsity pattern,
        // while others don't
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            var result = new double[Rows, Columns];
            for (int r = 0; r < Rows; ++r)
            {
                for (int c = 0; c < Columns; ++c)
                {
                    result[r, c] = unaryOperation(data[r, c]);
                }
            }
        }

        // Unsafe, but useful for passing the raw array to 3rd party libraries or doing optimizations after casting 
        // IMatrixView to DenseMatrix. I don't think you can mess the matrix more than you could by using the indexer.
        public double[,] ExposeInternalArray()
        {
            return data;
        }

        public IVector ExtractColumn(int col)
        {
            double[] colVector = new double[Rows];
            for (int r = 0; r < Rows; ++r) colVector[r] = data[r, col];
            return DenseVector.CreateFromArray(colVector, false);
        }

        public IVector ExtractDiagonal() //Only for diagonal
        {
            double[] diagonal = new double[Rows];
            for (int i = 0; i < Rows; ++i) diagonal[i] = data[i, i];
            return DenseVector.CreateFromArray(diagonal, false);
        }

        public IVector ExtractRow(int row)
        {
            double[] rowVector = new double[Columns];
            for (int c = 0; c < Columns; ++c) rowVector[c] = data[row, c];
            return DenseVector.CreateFromArray(rowVector, false);
        }

        public IMatrix ExtractSubmatrix(int[] rows, int[] columns)
        {
            var submatrix = new double[rows.Length, columns.Length];
            for (int r = 0; r < rows.Length; ++r)
            {
                int row = rows[r];
                for (int c = 0; c < columns.Length; ++c)
                {
                    submatrix[r, c] = data[row, columns[c]];
                }
            }
            return new DenseMatrix(submatrix);
        }

        public IVector MultiplyLeft(IVectorView vector)
        {
            Preconditions.CheckMultiplicationDimensions(vector, this);
            var result = new double[Columns];
            for (int c = 0; c < Columns; ++c) //Perhaps other access patterns are more efficient
            {
                double sum = 0.0;
                for (int r = 0; r < Rows; ++r)
                {
                    sum += vector[r] * data[r, c];
                }
                result[c] = sum;
            }
            return DenseVector.CreateFromArray(result, false);
        }

        public IVector MultiplyRight(IVectorView vector)
        {
            Preconditions.CheckMultiplicationDimensions(this, vector);
            var result = new double[Rows];
            for (int r = 0; r < Rows; ++r)
            {
                double sum = 0.0;
                for (int c = 0; c < Columns; ++c)
                {
                    sum += data[r, c] * vector[c];
                }
                result[r] = sum;
            }
            return DenseVector.CreateFromArray(result, false);
        }

        public IMatrix MultiplyRight(IMatrixView matrix)
        {
            Preconditions.CheckMultiplicationDimensions(this, matrix);
            var result = new double[this.Rows, matrix.Columns];
            for (int r = 0; r < this.Rows; ++r)
            {
                for (int c = 0; c < matrix.Columns; ++c)
                {
                    double sum = 0.0;
                    for (int k = 0; k < this.Columns; ++k)
                    {
                        sum += this.data[r, k] * matrix[k, c]; //perhaps if other is a DenseMatrix, matrix.data[k,c] for less redirection.
                    }
                    result[r, c] = sum;
                }
            }
            return new DenseMatrix(result);
        }

        public void Print()
        {
            for (int r = 0; r < Rows; ++r)
            {
                for (int c = 0; c < Columns; ++c)
                {
                    Console.Write(data[r, c]);
                    Console.Write(' '); // The excess space in last column doesn't hurt afaik.
                }
                Console.WriteLine();
            }
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, 
            ProcessZeros processZeros, Finalize finalize)
        {
            double accumulator = identityValue;
            for (int r = 0; r < Rows; ++r)
            {
                for (int c = 0; c < Columns; ++c)
                {
                    accumulator = processEntry(data[r, c], accumulator);
                }
            }
            // no zeros implied
            return finalize(accumulator);
        }

        public void SetAll(double value)
        {
            for (int r = 0; r < Rows; ++r)
            {
                for (int c = 0; c < Columns; ++c)
                {
                    data[r, c] = value;
                }
            }
        }

        public IMatrix Transpose()
        {
            var transpose = new double[Columns, Rows];
            for (int r = 0; r < Rows; ++r)
            {
                for (int c = 0; c < Columns; ++c)
                {
                    transpose[c, r] = data[r, c];
                }
            }
            return new DenseMatrix(transpose);
        }

        public IMatrixView TransposedView()
        {
            return new TransposedDenseMatrix(this, data);
        }

        public void Write(string path)
        {
            var writer = new StreamWriter(path);
            for (int r = 0; r < Rows; ++r)
            {
                for (int c = 0; c < Columns; ++c)
                {
                    writer.Write(data[r, c]);
                    writer.Write(' '); // The excess space in last column doesn't hurt afaik.
                }
                writer.WriteLine();
            }
#if DEBUG
            writer.Flush();
#endif
        }
    }
}
