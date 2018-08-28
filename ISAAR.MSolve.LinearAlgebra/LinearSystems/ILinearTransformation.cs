using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Not a big fan of the generics here. Perhaps I can avoid them by 1) having the caller pass in an instance of the result 
//      vector and writing into that, 2) Defining vector axpy and dot operations in an interface similar to this one. 
//      Alternatively I could break down the vector interfaces to smaller and more and then gather only the ones needed here.
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems
{
    /// <summary>
    /// Defines matrix-vector multiplication to allow iterative algorithms to operate without any modifications on various 
    /// matrix and vector types, such as distributed matrices and vectors.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TVector">The vector type that can be multiplied with this <see cref="ILinearTransformation{TVector}"/>.
    ///     </typeparam>
    public interface ILinearTransformation<TVector>
        where TVector: IVectorView
    {
        /// <summary>
        /// Performs the matrix-vector multiplication: this * <paramref name="vector"/>.
        /// </summary>
        /// <param name="vector">The vector that will be multiplied with this <see cref="ILinearTransformation{TVector}"/>.
        ///     </param>
        TVector Multiply(TVector vector); //TODO: should this be named Transform?
    }
}
