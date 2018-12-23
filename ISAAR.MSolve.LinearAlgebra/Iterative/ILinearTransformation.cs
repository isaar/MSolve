using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative
{
    /// <summary>
    /// Defines matrix-vector multiplication to allow iterative algorithms to operate without any modifications on various 
    /// matrix and vector types, such as distributed matrices and vectors.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ILinearTransformation
    {
        /// <summary>
        /// The number of columns of the matrix represented by this <see cref="ILinearTransformation"/>.
        /// </summary>
        int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix represented by this <see cref="ILinearTransformation"/>.
        /// </summary>
        int NumRows { get; }

        /// <summary>
        /// Performs the matrix-vector multiplication (with the matrix represented by this 
        /// <see cref="ILinearTransformation"/>): <paramref name="rhsVector"/> = this * <paramref name="lhsVector"/>.
        /// </summary>
        /// <param name="lhsVector">
        /// The vector that will be multiplied by the represented matrix. It sits on the left hand side of the equation 
        /// y = A * x. Constraints: Its <see cref="IIndexable1D.Length"/> must be equal to the number of columns of the matrix  
        /// represented by this <see cref="ILinearTransformation"/>.
        /// </param>
        /// <param name="rhsVector">
        /// The vector that will be overwritten by the result of the multiplication. It sits on the right hand side of the 
        /// equation y = A * x. Constraints: Its <see cref="IIndexable1D.Length"/> must be equal to the number of rows of the
        /// matrix represented by this <see cref="ILinearTransformation"/>.
        /// </param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="lhsVector"/> or <paramref name="rhsVector"/> violate the described constraints.
        /// </exception>
        void Multiply(IVectorView lhsVector, IVector rhsVector);
    }
}
