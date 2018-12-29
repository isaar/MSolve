namespace ISAAR.MSolve.LinearAlgebra.Iterative.Termination
{
    /// <summary>
    /// Determines the maximum number of iterations that will be run by an iterative algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IMaxIterationsProvider
    {
        /// <summary>
        /// Returns the max number of iterations appropriate for the specified <paramref name="matrix"/> of the linear system.
        /// </summary>
        /// <param name="matrix">The linear system's matrix. Must be square matrix.</param>
        int GetMaxIterationsForMatrix(ILinearTransformation matrix);
    }
}
