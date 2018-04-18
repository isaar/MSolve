using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// It supports common operations that do not mutate the underlying matrix. If you need to store a matrix and then pass it
    /// around or allow acceess to it, consider using this interface instead of <see cref="Matrix"/> for extra safety.
    /// </summary>
    public interface IMatrixView: IIndexable2D, IReducible
    {//TODO: Perhaps Addition, Subtraction and Scaling must be done without using delegates, for performance
        /// <summary>
        /// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. Optimized 
        /// version of <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/> and 
        /// <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>. Named after BLAS axpy (y = a*x plus y).
        /// </summary>
        /// <param name="other"></param>
        /// <param name="otherCoefficient"></param>
        /// <returns></returns>
        IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient);

        
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


        /// <summary>
        /// result[i, j] = <paramref name="thisCoefficient"/> * this[i, j] + <paramref name="otherCoefficient"/> * 
        /// <paramref name="otherMatrix"/>[i, j]. Optimized version of 
        /// <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        /// <param name="thisCoefficient"></param>
        /// <param name="other"></param>
        /// <param name="otherCoefficient"></param>
        /// <returns></returns>
        IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient);

        Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false);
        Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false);
        Vector MultiplyRight(IVectorView vector, bool transposeThis = false);
        IMatrixView Transpose();
    }
}
