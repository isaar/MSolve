using ISAAR.MSolve.LinearAlgebra.Vectors;

// Inversion is best handled by the matrix object itself, since the original should overwrite the factorized data in most cases, 
// which should be hidden from the user. Besides, I am not sure if first factorizing the matrix is more efficient than 
// Gauss-Jordan.
//TODO: also specify dimensions. For now this triangulation only works for square matrices.
namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// Represents the factorization of a square matrix into triangular and possible diagonal matrices. Such factorizations are 
    /// used to solve linear systems.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ITriangulation
    {
        /// <summary>
        /// Calculates the determinant of the original matrix (before the factorization). 
        /// </summary>
        double CalcDeterminant();

        /// <summary>
        /// Solves the linear system A * x = b, where A is the original matrix (before the factorization), 
        /// b = <paramref name="rhs"/> and x is the solution vector, which will overwrite the provided 
        /// <paramref name="solution"/>.
        /// </summary>
        /// <param name="rhs">
        /// The right hand side vector. Its <see cref="IIndexable1D.Length"/> must be equal to 
        /// <see cref="Matrices.IIndexable2D.NumRows"/> of the original matrix A.
        /// </param>
        /// <param name="solution">
        /// Output vector that will be overwritten with the solution of the linear system. Its <see cref="IIndexable1D.Length"/>  
        /// must be equal to <see cref="Matrices.IIndexable2D.NumColumns"/> of the original matrix A.
        /// </param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">
        /// Thrown if <paramref name="rhs"/> or <paramref name="solution"/> violate the described constraints.
        /// </exception>
        void SolveLinearSystem(Vector rhs, Vector solution);
    }
}
