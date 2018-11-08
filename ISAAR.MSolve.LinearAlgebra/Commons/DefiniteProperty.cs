//TODO: Perhaps the negative cases are also needed.
namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Property of a symmetric matrix. For more see: <see href="https://en.wikipedia.org/wiki/Positive-definite_matrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public enum DefiniteProperty
    {
        /// <summary>
        /// We are sure that the matrix is positive definite.
        /// </summary>
        PositiveDefinite,

        /// <summary>
        /// We are sure that the matrix is positive semi-definite (including positive definite).
        /// </summary>
        PositiveSemiDefinite,

        /// <summary>
        /// We are sure that the matrix is not positive definite, positive semi-definite, negative definite or negative 
        /// semi-definite.
        /// </summary>
        Indefinite,

        /// <summary>
        /// We do not yet know the definiteness of the matrix.
        /// </summary>
        Unknown
    }
}
