using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.LinearSystems
{
    /// <summary>
    /// An objects that defines a matrix-vector multiplication. It allows iterative algorithms to operate on various
    /// matrix and vector types, such as distributed matrices and vectors, without any modifications.
    /// </summary>
    /// <typeparam name="TVector"></typeparam>
    public interface ILinearTransformation<TVector>
        where TVector: IVectorView
    {
        /// <summary>
        /// y = T(x)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        TVector Multiply(TVector x); //TODO: should this be named Transform?
    }
}
