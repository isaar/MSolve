//TODO: make sure these match the pinvoke method arguments
namespace ISAAR.MSolve.LinearAlgebra.SuiteSparse
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

}
