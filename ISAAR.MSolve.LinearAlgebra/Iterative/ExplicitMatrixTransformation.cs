using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: each matrix-vector multiplication will internally cast the vectors, in order to perform efficiently. Could I avoid all
//      that casting by using generics
namespace ISAAR.MSolve.LinearAlgebra.Iterative
{
    /// <summary>
    /// Wrapper for a matrix class so that it can be used by iterative algorithms, which operate on 
    /// <see cref="ILinearTransformation"/>
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ExplicitMatrixTransformation : ILinearTransformation
    {
        private readonly IMatrixView matrix;

        /// <summary>
        /// Initializes a new instance of <see cref="ExplicitMatrixTransformation"/> that wraps the provided <paramref name="matrix"/>.
        /// </summary>
        /// <param name="matrix">The matrix that will be multiplied with vectors during the iterative algorithms.</param>
        public ExplicitMatrixTransformation(IMatrixView matrix) => this.matrix = matrix;

        /// <summary>
        /// See <see cref="ILinearTransformation.NumColumns"/>
        /// </summary>
        public int NumColumns => matrix.NumColumns;

        /// <summary>
        /// See <see cref="ILinearTransformation.NumRows"/>
        /// </summary>
        public int NumRows => matrix.NumRows;

        /// <summary>
        /// See <see cref="ILinearTransformation.Multiply(IVectorView, IVector)"/>
        /// </summary>
        public void Multiply(IVectorView lhsVector, IVector rhsVector) => matrix.MultiplyIntoResult(lhsVector, rhsVector, false);
    }
}
