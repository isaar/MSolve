using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    /// Operations specified by this interface modify the matrix. Therefore it is possible that they fail if they are used on 
    /// sparse matrices.
    /// </summary>
    public interface IMatrix: IMatrixView
    {
        /// <summary>
        /// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. Optimized 
        /// version of <see cref="IMatrix.DoEntrywise(IMatrixView, Func{double, double, double})"/> and 
        /// <see cref="IMatrix.LinearCombination(double, IMatrixView, double)"/>. Named after BLAS axpy (y = a*x plus y).
        /// </summary>
        /// <param name="other"></param>
        /// <param name="otherCoefficient"></param>
        /// <returns></returns>
        void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient);

        /// <summary>
        /// this[i, j] = binaryOperation(this[i,j], other[i,j])
        /// </summary>
        /// <param name="other">The other matrix. Must have the same dimensions as this instance.</param>
        /// <param name="binaryOperation">A function taking two arguments and returning one output.</param>
        /// <returns></returns>
        void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation);

        /// <summary>
        /// this[i, j] = unaryOperation(this[i,j])
        /// </summary>
        /// <param name="unaryOperation"></param>
        /// <returns></returns>
        void DoToAllEntriesIntoThis(Func<double, double> unaryOperation);

        /// <summary>
        /// this[i, j] = <paramref name="thisCoefficient"/> * this[i, j] + <paramref name="otherCoefficient"/> * 
        /// <paramref name="otherMatrix"/>[i, j]. Optimized version of 
        /// <see cref="IMatrix.DoEntrywise(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        /// <param name="thisCoefficient"></param>
        /// <param name="other"></param>
        /// <param name="otherCoefficient"></param>
        /// <returns></returns>
        void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient);

        /// <summary>
        /// Will throw a <see cref="Exceptions.SparsityPatternModifiedException"/> if a structural zero entry is written to.
        /// For symmetric matrices, this will set both (<paramref name="rowIdx"/>, <paramref name="colIdx"/>) and 
        /// (<paramref name="colIdx"/>, <paramref name="rowIdx"/>).
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        /// <param name="value"></param>
        void SetEntryRespectingPattern(int rowIdx, int colIdx, double value);
    }
}
