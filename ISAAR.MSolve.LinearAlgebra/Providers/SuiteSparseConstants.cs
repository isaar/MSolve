//TODO: make sure these match the pinvoke method arguments
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Settings to control the ordering of the matrix during factorization with the SuiteSparse library.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public enum SuiteSparseOrdering
    {
        /// <summary>
        /// No reordering. This allows the most control, as you can first find a fill-reducing ordering for the sparsity pattern,
        /// apply it and finally factorize using <see cref="SuiteSparseOrdering.Natural"/>.
        /// </summary>
        Natural = 0,

        /// <summary>
        /// Let SuiteSparse try all its default algorithms, until it finds a good ordering.
        /// </summary>
        TryDefaults = 1,

        /// <summary>
        /// Use Approximate Minimal Degree algorithm.
        /// </summary>
        AMD = 2
    }

    /// <summary>
    /// How to proceed after factorizing (triangulation) the linear system matrix: A=L*U.
    /// </summary>
    internal enum SystemType
    {
        /// <summary>
        /// Solve the original system A*x = b
        /// </summary>
        Regular = 0,

        /// <summary>
        /// Solve L*x = b
        /// </summary>
        ForwardSubstitution = 4,

        /// <summary>
        /// Solve U*x = b
        /// </summary>
        BackSubstitution = 5
    }
}
