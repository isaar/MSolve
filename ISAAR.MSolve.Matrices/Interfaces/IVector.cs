using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Matrices.Interfaces
{
    public interface IVector<T>
    {
        int Length { get; }
        double Norm { get; }
        T this[int x] { get; set; }
        void Multiply(double coefficient);
        void CopyTo(Array array, int index);
        void Clear();
        void WriteToFile(string name);
    }
}
