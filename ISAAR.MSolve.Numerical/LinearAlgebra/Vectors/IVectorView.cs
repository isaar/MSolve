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
    /// <summary>
    /// It supports common operations that do not mutate the underlying vector. If you need to store a vector and then pass it
    /// around or allow acceess to it, consider using this interface instead of <see cref="VectorMKL"/> for extra safety.
    /// </summary>
    public interface IVectorView: IReducible
    {
        double this[int index] { get; }
        int Length { get; }
      
        double[] CopyToArray();
        IVectorView DoEntrywise(IVectorView other, Func<double, double, double> binaryOperation);
        IVectorView DoToAllEntries(Func<double, double> unaryOperation);
        double DotProduct(IVectorView vector);
    }
}
