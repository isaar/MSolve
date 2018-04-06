using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    /// It supports common operations that do not mutate the underlying matrix. If you need to store a matrix and then pass it
    /// around or allow acceess to it, consider using this interface instead of <see cref="Matrix"/> for extra safety.
    /// </summary>
    public interface IMatrixView: IIndexable2D, IReducible
    {
        //TODO: Perhaps Addition, Subtraction and Scaling must be done without using delegates, for performance
        /// <summary>
        /// result[i, j] = binaryOperation(this[i,j], other[i,j])
        /// </summary>
        /// <param name="other">The other matrix. Must have the same dimensions as this instance.</param>
        /// <param name="binaryOperation">A function taking two arguments and returning one output.</param>
        /// <returns></returns>
        IMatrixView DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation);

        /// <summary>
        /// result[i, j] = unaryOperation(this[i,j])
        /// </summary>
        /// <param name="unaryOperation"></param>
        /// <returns></returns>
        IMatrixView DoToAllEntries(Func<double, double> unaryOperation);


        Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false);
        Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false);
        VectorMKL MultiplyRight(IVectorView vector, bool transposeThis = false);
        IMatrixView Transpose();
    }
}
