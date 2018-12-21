using IntelMKL.LP64;

//TODO: this should probably call the MKL dll directly, instead of using the package Compute.NET Bindings.
namespace ISAAR.MSolve.LinearAlgebra.Providers.MKL
{
    internal class MklBlasExtensionsProvider : IBlasExtensionsProvider
    {
        internal static MklBlasExtensionsProvider UniqueInstance { get; } = new MklBlasExtensionsProvider();

        private MklBlasExtensionsProvider() { } // private constructor for singleton pattern

        public void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY,
            int incY)
            => Blas.Daxpby(ref n, ref alpha, ref x[offsetX], ref incX, ref beta, ref y[offsetY], ref incY);
    }
}
