using System.Collections.Generic;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces
{
    public interface IMatrix2D
    {
        int Rows { get; }
        int Columns { get; }
        double this[int x, int y] { get; set; }
        void Scale(double scale);
        void Multiply(IVector vIn, double[] vOut);
    }
}
