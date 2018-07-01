using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.IO;
using System.Globalization;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class SymmetricMatrix2D : IMatrix2D, ILinearlyCombinable, ILinearlyCombinable<SymmetricMatrix2D>
    {
        private int rows;
        private double[] data;

        public SymmetricMatrix2D(int rows)
        {
            this.rows = rows;
            data = new double[(rows + 1) * rows / 2];
        }

        public SymmetricMatrix2D(double[] data)
        {
            this.rows = (int)(Math.Sqrt(8 * data.Length + 1) - 1) / 2;
            this.data = data;
        }

        public SymmetricMatrix2D(Matrix2D matrix)
        {
            if (matrix.Columns != matrix.Rows)
                throw new ArgumentException("Matrix2D is NOT rectangular.");
            
            this.rows = matrix.Rows;
            data = new double[(rows + 1) * rows / 2];
            int count = 0;
            for (int i = 0; i < rows; i++)
                for (int j = i; j < rows; j++)
                {
                    data[count] = matrix[i, j];
                    count++;
                }
        }

        public int IndexOf(int row, int column)
        {
            int r = Math.Min(row, column);
            int c = Math.Max(row, column);
            return data.Length - ((rows - r + 1) * (rows - r) / 2) + c - r;
        }

        public double[] Data
        {
            get { return data; }
        }

        #region IMatrix2D<T> Members

        public int Rows
        {
            get { return rows; }
        }

        public int Columns
        {
            get { return rows; }
        }

        public double this[int x, int y]
        {
            get { return data[IndexOf(x, y)]; }
            set { data[IndexOf(x, y)] = value; }
        }

        public static double[] operator *(SymmetricMatrix2D K, IVector v)
        {
            if (K.Rows != v.Length) throw new ArgumentException("Matrix and vector size mismatch.");
            double[] result = new double[K.Rows];
            K.Multiply(v, result);
            return result;
        }

        public Matrix2D ToMatrix2D()
        {
            var matrix = new double[rows, rows];
            int index = 0;
            for (int i = 0; i < rows; i++)
                for (int j = i; j < rows; j++)
                {
                    matrix[i, j] = data[index];
                    matrix[j, i] = data[index];
                    index++;
                }

            return new Matrix2D(matrix);
        }

        public void Multiply(IVector vIn, double[] vOut)
        {
            if (Rows != vIn.Length) throw new ArgumentException("Matrix and vector size mismatch.");

            Array.Clear(vOut, 0, vOut.Length);
            double[] d = data;
            //int pos = 0;
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                    vOut[i] += d[IndexOf(i, j)] * vIn[j];
        }

        public void Scale(double scale)
        {
            double[] mData = data;
            for (int i = 0; i < mData.Length; i++) mData[i] *= scale;
        }

        public void WriteToFile(string name)
        {
            double[] mData = data;

            StreamWriter sw = new StreamWriter(name);
            foreach (double d in mData)
                sw.WriteLine(d.ToString("g17", new CultureInfo("en-US", false).NumberFormat));
            sw.Close();
        }

        public void ReadFromFile(string name)
        {
            string[] lines = File.ReadAllLines(name);
            data = new double[lines.Length];
            double[] mData = data;
            for (int i = 0; i < lines.Length; i++) mData[i] = Convert.ToDouble(lines[i]);
        }

        public void LinearCombination(IList<double> coefficients, IList<SymmetricMatrix2D> matrices)
        {
            if (coefficients.Count != matrices.Count)
                throw new ArgumentException(String.Format("Coefficients and matrices count mismatch ({0} <> {1}).", coefficients.Count, matrices.Count));
            for (int i = 0; i < matrices.Count; i++)
                if (matrices[i].Rows != rows)
                    throw new ArgumentException(String.Format("Matrix at pos {0} has {1} rows instead of {2}.", i, matrices[i].Rows, rows));

            var cs = (IList<double>)coefficients;
            double[] newData = new double[data.Length];
            for (int i = 0; i < matrices.Count; i++)
            {
                var m = matrices[i];
                for (int j = 0; j < data.Length; j++)
                    newData[j] += cs[i] * m.Data[j];
            }

            var d = data;
            Array.Copy(newData, d, data.Length);
        }

        public void LinearCombination(IList<double> coefficients, IList<IMatrix2D> matrices)
        {
            if (coefficients.Count != matrices.Count)
                throw new ArgumentException(String.Format("Coefficients and matrices count mismatch ({0} <> {1}).", coefficients.Count, matrices.Count));
            for (int i = 0; i < matrices.Count; i++)
                if (matrices[i].Rows != rows)
                    throw new ArgumentException(String.Format("Matrix at pos {0} has {1} rows instead of {2}.", i, matrices[i].Rows, rows));

            var cs = coefficients.ToArray();
            double[] newData = new double[data.Length];
            for (int k = 0; k < matrices.Count; k++)
            {
                var m = matrices[k];
                int index = 0;
                for (int i = 0; i < rows; i++)
                    for (int j = i; j < rows; j++)
                        newData[index++] += cs[k] * m[i,j];
                //for (int j = 0; j < data.Length; j++)
                //    newData[j] += cs[i] * m.Data[j];
            }

            var d = data;
            Array.Copy(newData, d, data.Length);
        }

        #endregion
    }
}
