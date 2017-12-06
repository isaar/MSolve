using System;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces
{
    public interface IVectorOLD
    {
        int Length { get; }
        double Norm { get; }
        double this[int x] { get; set; }
        double DotProduct(IVectorOLD y);
        void Multiply(double coefficient);
        void CopyTo(Array array, int index);
        void CopyFrom(int startIndex, int length, IVectorOLD fromVector, int fromStartIndex);
        void Clear();
    }
}
