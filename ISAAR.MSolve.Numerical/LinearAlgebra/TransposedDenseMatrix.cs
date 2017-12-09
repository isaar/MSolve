using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    //Perhaps this could be a private DenseMatrix class and delegate most of its operations there
    class TransposedDenseMatrix: IMatrixView
    {
        private readonly double[,] data; //shared with another matrix instance
        private readonly DenseMatrix untransposedMatrix;

        public double this[int row, int col]
        {
            get { return data[col, row]; }
            set { data[col, row] = value; }
        }

        public int Rows { get; }
        public int Columns { get; }
        public int NumNonZeros { get { return Rows * Columns; } }

        public TransposedDenseMatrix(DenseMatrix untransposedMatrix, double[,] dataOfUntransposedMatrix)
        {
            this.data = dataOfUntransposedMatrix;
            this.untransposedMatrix = untransposedMatrix;
            Rows = untransposedMatrix.Columns;
            Columns = untransposedMatrix.Rows;
        }

        public IMatrix DoPointwise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckMatrixDimensionsSame(this, other);
            var result = new double[this.Rows, this.Columns]; // These are correct. Only data is transposed.
            for (int r = 0; r < this.Rows; ++r)
            {
                for (int c = 0; c < this.Columns; ++c)
                {
                    result[r, c] = binaryOperation(this.data[c, r], other[r, c]);
                }
            }
            return DenseMatrix.CreateFromArray(result, false);
        }


        public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
        {
            return untransposedMatrix.DoToAllEntries(unaryOperation);
        }

        public IVector ExtractColumn(int col)
        {
            return untransposedMatrix.ExtractRow(col);
        }

        public IVector ExtractDiagonal()
        {
            return untransposedMatrix.ExtractDiagonal();
        }

        public IVector ExtractRow(int row)
        {
            return untransposedMatrix.ExtractColumn(row);
        }

        public IMatrix ExtractSubmatrix(int[] rows, int[] columns)
        {
            var submatrix = new double[rows.Length, columns.Length];
            for (int r = 0; r < rows.Length; ++r) // Perhaps iterating columns first is more cache-friendly
            {
                int row = rows[r];
                for (int c = 0; c < columns.Length; ++c)
                {
                    submatrix[r, c] = data[columns[c], row];
                }
            }
            return DenseMatrix.CreateFromArray(submatrix, false);
        }

        public IVector MultiplyLeft(IVectorView vector)
        {
            return untransposedMatrix.MultiplyRight(vector); // x^T*A^T = (A*x)^T
        }

        public IVector MultiplyRight(IVectorView vector)
        {
            return untransposedMatrix.MultiplyLeft(vector); // A^T*x = (x^T*A)^T
        }
        
        // Column major iteration of this matrix == row major iteration of this.data. Should this be in the IMatrixView?
        public DenseMatrix MultiplyLeft(IMatrixView matrix)
        {
            Preconditions.CheckMultiplicationDimensions(matrix, this);
            var result = new double[this.Rows, matrix.Columns];
            for (int r = 0; r < matrix.Rows; ++r)
            {
                for (int c = 0; c < this.Columns; ++c)
                {
                    double sum = 0.0;
                    for (int k = 0; k < matrix.Columns; ++k)
                    {
                        sum += matrix[r, k] * this.data[c, k]; //perhaps if other is a DenseMatrix, matrix.data[k,c] for less redirection.
                    }
                    result[r, c] = sum;
                }
            }
            return DenseMatrix.CreateFromArray(result, false);
        }

        // This will have row major accesses on this, thus column major accesses on this.data and column major accesses
        // on matrix. Performace is expected to be low.
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
                        sum += this.data[k, r] * matrix[k, c]; //perhaps if other is a DenseMatrix, matrix.data[k,c] for less redirection.
                    }
                    result[r, c] = sum;
                }
            }
            return DenseMatrix.CreateFromArray(result, false);
        }

        public void Print()
        {
            for (int r = 0; r < Rows; ++r)
            {
                for (int c = 0; c < Columns; ++c)
                {
                    Console.Write(data[c, r]);
                    Console.Write(' '); // The excess space in last column doesn't hurt afaik.
                }
                Console.WriteLine();
            }
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, 
            ProcessZeros processZeros, Finalize finalize)
        {
            // Order of entries does not matter.
            return untransposedMatrix.Reduce(identityValue, processEntry, processZeros, finalize); 
        }

        public IMatrix Transpose()
        {
            return DenseMatrix.CreateFromArray(data, true);
        }

        public IMatrixView TransposedView()
        {
            return DenseMatrix.CreateFromArray(data, false);
        }

        public void Write(string path)
        {
            var writer = new StreamWriter(path);
            for (int r = 0; r < Rows; ++r)
            {
                for (int c = 0; c < Columns; ++c)
                {
                    writer.Write(data[c, r]);
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
