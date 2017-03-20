using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;
using System.IO;
using System.Globalization;

namespace ISAAR.MSolve.Matrices
{
    public class SymmetricMatrix2D<T> : IMatrix2D<T>
    {
        private int rows;
        private T[] data;

        public SymmetricMatrix2D(int rows)
        {
            this.rows = rows;
            data = new T[(rows + 1) * rows / 2];
        }

        public SymmetricMatrix2D(T[] data)
        {
            this.rows = (int)(Math.Sqrt(8 * data.Length + 1) - 1) / 2;
            this.data = data;
        }

        public SymmetricMatrix2D(Matrix2D<T> matrix)
        {
            if (matrix.Columns != matrix.Rows)
                throw new ArgumentException("Matrix2D is NOT rectangular.");

            this.rows = matrix.Rows;
            data = new T[(rows + 1) * rows / 2];
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

        public T[] Data
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

        public T this[int x, int y]
        {
            get { return data[IndexOf(x, y)]; }
            set { data[IndexOf(x, y)] = value; }
        }

        public static double[] operator *(SymmetricMatrix2D<T> K, IVector<double> v)
        {
            if (!(typeof(T) == typeof(double))) throw new InvalidOperationException("Cannot solve for types other than double");
            if (K.Rows != v.Length) throw new InvalidOperationException("Matrix and vector size mismatch.");
            double[] result = new double[K.Rows];
            K.Multiply(v, result);
            return result;
        }

        public Matrix2D<T> ToMatrix2D()
        {
            T[,] matrix = new T[rows, rows];
            int index = 0;
            for (int i = 0; i < rows; i++)
                for (int j = i; j < rows; j++)
                {
                    matrix[i, j] = data[index];
                    matrix[j, i] = data[index];
                    index++;
                }

            return new Matrix2D<T>(matrix);
        }

        public void Solve(IVector<double> f, double[] result)
        {
            throw new NotImplementedException();
        }

        public void LinearCombination(IList<T> coefficients, IList<IMatrix2D<T>> matrices)
        {
            if (!(typeof(T) == typeof(double))) throw new InvalidOperationException("Cannot solve for types other than double");
            if (coefficients.Count != matrices.Count)
                throw new ArgumentException(String.Format("Coefficients and matrices count mismatch ({0} <> {1}).", coefficients.Count, matrices.Count));
            for (int i = 0; i < matrices.Count; i++)
            {
                if (matrices[i] is SymmetricMatrix2D<T> == false)
                    throw new ArgumentException(String.Format("Matrix at pos {0} is of type {1} instead of SymmetricMatrix2D.", i, matrices[i].GetType()));
                if (matrices[i].Rows != rows)
                    throw new ArgumentException(String.Format("Matrix at pos {0} has {1} rows instead of {2}.", i, matrices[i].Rows, rows));
            }

            var cs = (IList<double>)coefficients;
            double[] newData = new double[data.Length];
            for (int i = 0; i < matrices.Count; i++)
            {
                var m = matrices[i] as SymmetricMatrix2D<double>;
                for (int j = 0; j < data.Length; j++)
                    newData[j] += cs[i] * m.Data[j];
            }

            var d = data as double[];
            Array.Copy(newData, d, data.Length);
        }

        public void Multiply(IVector<double> vIn, double[] vOut)
        {
            if (!(typeof(T) == typeof(double))) throw new InvalidOperationException("Cannot multiply for types other than double");
            if (Rows != vIn.Length) throw new InvalidOperationException("Matrix and vector size mismatch.");

            Array.Clear(vOut, 0, vOut.Length);
            double[] d = data as double[];
            //int pos = 0;
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                    vOut[i] += d[IndexOf(i, j)] * vIn[j];
        }

        public void Scale(double scale)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] mData = data as double[];
            for (int i = 0; i < mData.Length; i++) mData[i] *= scale;
        }

        public void WriteToFile(string name)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] mData = data as double[];

            StreamWriter sw = new StreamWriter(name);
            foreach (double d in mData)
                sw.WriteLine(d.ToString("g17", new CultureInfo("en-US", false).NumberFormat));
            sw.Close();
        }

        public void ReadFromFile(string name)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");

            string[] lines = File.ReadAllLines(name);
            data = new T[lines.Length];
            double[] mData = data as double[];
            for (int i = 0; i < lines.Length; i++) mData[i] = Convert.ToDouble(lines[i]);
        }

        #endregion

        public override String ToString()
        {
            StringBuilder builder = new StringBuilder();
            for (int row = 0; row < Rows; ++row)
            {
                for (int col = 0; col < Columns; ++col)
                {
                    builder.Append(this[row, col]);
                    builder.Append(' ');
                }
                builder.Append("\n");
            }
            builder.Append("\n");
            return builder.ToString();
        }

        public String UpperTriangleToString()
        {
            StringBuilder builder = new StringBuilder();
            for (int row = 0; row < Rows; ++row)
            {
                for (int col = 0; col < Columns; ++col)
                {
                    if (col < row) builder.Append(0.0);
                    else builder.Append(this[row, col]);
                    builder.Append(' ');
                }
                builder.Append("\n");
            }
            builder.Append("\n");
            return builder.ToString();
        }
    }
}
