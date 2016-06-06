using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Matrices
{
    public class Sparse2D<T> : IMatrix2D<T>
    {
        private int rows, columns;
        private readonly Dictionary<int, Dictionary<int, int>> rowColDataPositions;
        private readonly List<T> data;

        public Sparse2D(int rows, int columns)
        {
            this.rows = rows;
            this.columns = columns;
            rowColDataPositions = new Dictionary<int, Dictionary<int, int>>(rows);
            data = new List<T>(rows);
        }

        private T GetValueFromRowCol(int row, int col)
        {
            T value = default(T);
            if (rowColDataPositions.ContainsKey(row))
                if (rowColDataPositions[row].ContainsKey(col)) value = data[rowColDataPositions[row][col]];
            return value;
        }

        private void SetValueAtRowCol(int row, int col, T value)
        {
            int pos = data.Count;
            if (rowColDataPositions.ContainsKey(row))
            {
                if (rowColDataPositions[row].ContainsKey(col))
                {
                    data[rowColDataPositions[row][col]] = value;
                    return;
                }
                else
                    rowColDataPositions[row].Add(col, pos);
            }
            else
            {
                rowColDataPositions.Add(row, new Dictionary<int, int>());
                rowColDataPositions[row].Add(col, pos);
            }
            data.Add(value);
        }

        #region IMatrix2D<double> Members

        public int Rows
        {
            get { return rows; }
        }

        public int Columns
        {
            get { return columns; }
        }

        public T this[int x, int y]
        {
            get { return GetValueFromRowCol(x, y); }
            set { SetValueAtRowCol(x, y, value); }
        }

        public void Scale(double scale)
        {
            throw new NotImplementedException();
        }

        public void Multiply(IVector<double> vIn, double[] vOut)
        {
            if (!(typeof(T) == typeof(double))) throw new InvalidOperationException("Cannot multiply for types other than double");
            if (Rows != vIn.Length) throw new InvalidOperationException("Matrix and vector size mismatch.");
            Array.Clear(vOut, 0, vOut.Length);
            List<double> d = data as List<double>;

            foreach (int row in rowColDataPositions.Keys)
                foreach (int col in rowColDataPositions[row].Keys)
                    vOut[row] += d[rowColDataPositions[row][col]] * vIn[col];
        }

        public void LinearCombination(IList<T> coefficients, IList<IMatrix2D<T>> matrices)
        {
            throw new NotImplementedException();
        }

        public void Solve(IVector<double> f, double[] result)
        {
            throw new NotImplementedException();
        }

        public void WriteToFile(string name)
        {
            throw new NotImplementedException();
        }

        public void ReadFromFile(string name)
        {
            throw new NotImplementedException();
        }

        #endregion
    }
}
