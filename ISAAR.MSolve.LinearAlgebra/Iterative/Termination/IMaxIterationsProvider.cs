namespace ISAAR.MSolve.LinearAlgebra.Iterative.Termination
{
    /// <summary>
    /// Determines the maximum number of iterations that will be run by an iterative algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IMaxIterationsProvider
    {
        /// <summary>
        /// Returns the max number of iterations appropriate for the specified order of the linear system's matrix.
        /// </summary>
        /// <param name="matrix">The order of the linear system's matrix.</param>
        int GetMaxIterations(int matrixOrder);
    }
}
