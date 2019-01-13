//TODO: make sure these match the pinvoke method arguments
//TODO: These reorderings must be available to the user as dedicated classes. After selecting one the dofs of the finite 
//      element model should be reordered, since that reordering will be used many times. CHOLDMOD should always be used with 
//      "natural" ordering, which may have been produced by the user by applying a reordering algorithm.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Settings to control the ordering of the matrix during factorization with the SuiteSparse library.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal enum SuiteSparseOrdering 
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
    /// Specifies whether to solve the linear system A * x = b with the original matrix A or with some factor produced by the 
    /// factorization A = L * U or A = L * transpose(L).
    /// Authors: Serafeim Bakalakos
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
