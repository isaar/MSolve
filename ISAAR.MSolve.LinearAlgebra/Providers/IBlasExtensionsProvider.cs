namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Provides linear algebra operations that are variations of the usual BLAS subroutines (see <see cref="IBlasProvider"/>) 
    /// or in the same spirit.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal interface IBlasExtensionsProvider
    {
        /// <summary>
        /// y = alpha * x + beta * y
        /// </summary>
        void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY);
    }
}
