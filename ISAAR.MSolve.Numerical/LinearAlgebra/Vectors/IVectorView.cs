using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Vectors
{
    public interface IVectorView: IReducible
    {
        double this[int index] { get; }
        int Length { get; }
      
        double[] CopyToArray();
        VectorMKL DoPointwise(IVectorView other, Func<double, double, double> operation);
        VectorMKL DoToAllEntries(Func<double, double> operation);
        double DotProduct(IVectorView vector);
    }
}
