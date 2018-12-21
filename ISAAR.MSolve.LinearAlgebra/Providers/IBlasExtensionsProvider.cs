namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    internal interface IBlasExtensionsProvider
    {
        /// <summary>
        /// y = alpha * x + beta * y
        /// </summary>
        void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY);
    }
}
