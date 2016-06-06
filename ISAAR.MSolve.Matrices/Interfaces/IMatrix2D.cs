using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Matrices.Interfaces
{
    public interface IMatrix2D<T>
    {
        int Rows { get; }
        int Columns { get; }
        T this[int x, int y] { get; set; }
        void Scale(double scale);
        void Multiply(IVector<double> vIn, double[] vOut);
        void LinearCombination(IList<T> coefficients, IList<IMatrix2D<T>> matrices);
        void Solve(IVector<double> f, double[] result);
        void WriteToFile(string name);
        void ReadFromFile(string name);
    }
}
