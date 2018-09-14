using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.LinearSystems
{
    /// <summary>
    /// Wrapper for a matrix class, so that it can be used by iterative algorithms, which operate on 
    /// <see cref="ILinearTransformation{TVector}"/>
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class MatrixTransformation: ILinearTransformation<Vector>
    {
        private readonly IMatrixView matrix;

        /// <summary>
        /// Initializes a new instance of <see cref="MatrixTransformation"/> that wraps the provided <paramref name="matrix"/>.
        /// </summary>
        /// <param name="matrix">The matrix that will be multiplied with vectors during the iterative algorithms.</param>
        public MatrixTransformation(IMatrixView matrix) => this.matrix = matrix;

        /// <summary>
        /// See <see cref="ILinearTransformation{TVector}.Multiply(TVector)"/>
        /// </summary>
        public Vector Multiply(Vector vector) => matrix.MultiplyRight(vector, false);
    }
}
