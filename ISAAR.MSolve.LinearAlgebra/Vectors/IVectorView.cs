using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// It supports common operations that do not mutate the underlying vector. If you need to store a vector and then pass it
    /// around or allow acceess to it, consider using this interface instead of <see cref="Vector"/> for extra safety.
    /// </summary>
    public interface IVectorView: IIndexable1D, IReducible
    {
        double[] CopyToArray();
        IVectorView DoEntrywise(IVectorView other, Func<double, double, double> binaryOperation);
        IVectorView DoToAllEntries(Func<double, double> unaryOperation);
        double DotProduct(IVectorView vector);
    }
}
