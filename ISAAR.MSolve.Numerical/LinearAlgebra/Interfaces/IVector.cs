using System;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces
{
    public interface IVector
    {
        int Length { get; }
        double Norm { get; }
        double this[int x] { get; set; }
        double DotProduct(IVector y);
        void Multiply(double coefficient);
        void CopyTo(Array array, int index);
        void CopyFrom(int startIndex, int length, IVector fromVector, int fromStartIndex);
        void Clear();
        Vector[] RemoveDuplicatesFindMultiplicity();
    }
}
