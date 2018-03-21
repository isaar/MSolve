using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public interface IVectorView: IReducible
    {
        double this[int index] { get; }
        int Length { get; }
      
        double[] CopyToArray();
        Vector DoPointwise(IVectorView other, Func<double, double, double> operation);
        Vector DoToAllEntries(Func<double, double> operation);
        double DotProduct(IVectorView vector);
        Vector Slice(int[] indices);
        Vector Slice(int startInclusive, int endExclusive);
        void WriteToConsole(Array1DFormatting format = null);
        void WriteToFile(string path, bool append = false, Array1DFormatting format = null);
    }
}
