using ISAAR.MSolve.LinearAlgebra.Vectors;

// Inversion is best handled by the matrix object itself, since the original should overwrite the factorized data in most cases, 
// which should be hidden from the user. Besides, I am not sure if first factorizing the matrix is more efficient than 
// Gauss-Jordan.
//TODO: also specify dimensions. For now this triangulation only works for square matrices.
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
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
        /// b = <paramref name="rhsVector"/> and x is the result of this method.
        /// </summary>
        /// <param name="rhsVector">A vector with <see cref="IIndexable1D.Length"/> equal to the 
        ///     <see cref="Matrices.IIndexable2D.NumRows"/> of the original matrix A.</param>
        Vector SolveLinearSystem(Vector rhsVector);
    }
}
