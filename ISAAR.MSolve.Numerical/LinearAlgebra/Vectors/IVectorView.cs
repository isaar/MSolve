using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Vectors
{
    public interface IVectorView: IReducible, IWriteable
    {
        double this[int index] { get; }
        int Length { get; }
      
        double[] CopyToArray();
        VectorMKL DoPointwise(IVectorView other, Func<double, double, double> operation);
        VectorMKL DoToAllEntries(Func<double, double> operation);
        double DotProduct(IVectorView vector);
        VectorMKL Slice(int[] indices);
        VectorMKL Slice(int startInclusive, int endExclusive);
    }
}
