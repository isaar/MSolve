using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class Sparse2D : IMatrix2D
    {
        private int rows, columns;
        private readonly Dictionary<int, Dictionary<int, int>> rowColDataPositions;
        private readonly List<double> data;

        public Sparse2D(int rows, int columns)
        {
            this.rows = rows;
            this.columns = columns;
            rowColDataPositions = new Dictionary<int, Dictionary<int, int>>(rows);
            data = new List<double>(rows);
        }

        private double GetValueFromRowCol(int row, int col)
        {
            return rowColDataPositions.ContainsKey(row) && rowColDataPositions[row].ContainsKey(col) ? data[rowColDataPositions[row][col]] : 0;
            //double value = 0;
            //if (rowColDataPositions.ContainsKey(row))
            //    if (rowColDataPositions[row].ContainsKey(col)) value = data[rowColDataPositions[row][col]];
            //return value;
        }

        private void SetValueAtRowCol(int row, int col, double value)
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

        public double this[int x, int y]
        {
            get { return GetValueFromRowCol(x, y); }
            set { SetValueAtRowCol(x, y, value); }
        }

        public void Scale(double scale)
        {
            throw new NotImplementedException();
        }

        public void Multiply(IVector vIn, double[] vOut)
        {
            if (Columns != vIn.Length) throw new ArgumentException("Matrix and vector size mismatch.");
            Array.Clear(vOut, 0, vOut.Length);
            List<double> d = data as List<double>;

            foreach (int row in rowColDataPositions.Keys)
                foreach (int col in rowColDataPositions[row].Keys)
                    vOut[row] += d[rowColDataPositions[row][col]] * vIn[col];
        }

        #endregion
    }
}
