using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.LinearAlgebra
{
    /// <summary>
    /// A simple unoptimized immutable double[,] with some useful methods for linear algebra operations. 
    /// Use it for small matrices.
    /// </summary>
    class DenseMatrix
    {
        public static DenseMatrix CreateByCopying(double[,] data)
        {
            double[,] newData = new double[data.GetLength(0), data.GetLength(1)];
            Array.Copy(data, newData, data.Length);
            return new DenseMatrix(newData);
        }

        public static double[] operator *(DenseMatrix matrix, double[] columnVector)
        {
            return matrix.MultiplyRight(columnVector);
        }

        public static double[] operator *(double[] rowVector, DenseMatrix matrix)
        {
            return matrix.MultiplyLeft(rowVector);
        }

        public static DenseMatrix operator *(DenseMatrix leftMatrix, DenseMatrix rightMatrix)
        {
            return leftMatrix.MultiplyRight(rightMatrix);
        }

        private readonly int rows, columns;
        private readonly double[,] data;

        public double this[int row, int col]
        {
            get { return data[row, col]; }
        }

        public DenseMatrix(double[,] data)
        {
            this.rows = data.GetLength(0);
            this.columns = data.GetLength(1);
            this.data = data;
        }

        public DenseMatrix Transpose()
        {
            double[,] transpose = new double[rows, columns];
            for (int r = 0; r < rows; ++r)
            {
                for (int c = 0; c < columns; ++c)
                {
                    transpose[r, c] = data[c, r];
                }
            }
            return new DenseMatrix(transpose);
        }

        /// <summary>
        /// Performs the matrix-vector multiplication: this * columnVector. Doesn't check dimensions.
        /// </summary>
        /// <param name="columnVector">An array with length == the columns count of this matrix</param>
        /// <returns></returns>
        public double[] MultiplyRight(double[] columnVector)
        {
            CheckDimensionsThisTimesVector(columnVector);
            double[] result = new double[rows];
            for (int r = 0; r < rows; ++r)
            {
                double sum = 0.0;
                for (int c = 0; c < columns; ++c) sum += data[r, c] * columnVector[c];
                result[r] = sum;
            }
            return result;
        }

        /// <summary>
        /// Performs the following matrix multiplication: this * matrix. Doesn't check dimensions.
        /// </summary>
        /// <param name="otherMatrix">A matrix with rows count == the columns count of this matrix</param>
        /// <returns></returns>
        public DenseMatrix MultiplyRight(DenseMatrix otherMatrix)
        {
            CheckDimensionsThisTimesOtherMatrix(otherMatrix);
            double[,] result = new double[this.rows, otherMatrix.columns];
            for (int r = 0; r < this.rows; ++r)
            {
                for (int c = 0; c < otherMatrix.columns; ++c)
                {
                    double sum = 0.0;
                    for (int i = 0; i < this.columns; ++i) sum += this.data[r, i] * otherMatrix.data[i, c];
                    result[r, c] = sum;
                }
            }
            return new DenseMatrix(result);
        }

        /// <summary>
        /// Performs the following matrix-vector multiplication: rowVector * this. Doesn't check dimensions.
        /// </summary>
        /// <param name="matrix">An array with length == the rows count of this matrix</param>
        /// <returns></returns>
        public double[] MultiplyLeft(double[] rowVector)
        {
            CheckDimensionsVectorTimesThis(rowVector);
            double[] result = new double[columns];
            for (int c = 0; c < columns; ++c)
            {
                double sum = 0.0;
                for (int r = 0; r < rows; ++r) sum += data[r, c] * rowVector[r];
                result[c] = sum;
            }
            return result;
        }

        [System.Diagnostics.Conditional("DEBUG")]
        private void CheckDimensionsThisTimesVector(double[] columnVector)
        {
            if (columnVector.Length != columns) throw new NonMatchingDimensionsException("Matrix (" + rows 
                + ", " + columns + ") x Vector (" + columnVector.Length + ")");
        }

        [System.Diagnostics.Conditional("DEBUG")]
        private void CheckDimensionsVectorTimesThis(double[] rowVector)
        {
            if (rowVector.Length != rows) throw new NonMatchingDimensionsException("Vector (" + rowVector.Length +
                ") x Matrix (" + rows  + ", " + columns + ")");
        }

        [System.Diagnostics.Conditional("DEBUG")]
        private void CheckDimensionsThisTimesOtherMatrix(DenseMatrix other)
        {
            if (this.columns != other.rows) throw new NonMatchingDimensionsException("Matrix (" + this.rows
                + ", " + this.columns + ") x Matrix (" + other.rows + ", " + other.columns + ")");
        }
    }
}
