using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public interface IVectorView: IReducible
    {
        double this[int index] { get; }
        int Length { get; }
      
        double[] CopyToArray();
        IVector DoPointwise(IVectorView other, Func<double, double, double> operation);
        IVector DoToAllEntries(Func<double, double> operation);
        IVector ExtractSubvector(int[] indices);
        double MultiplyDot(IVectorView vector);
        void Print();
        void Write(string path);
    }
}
