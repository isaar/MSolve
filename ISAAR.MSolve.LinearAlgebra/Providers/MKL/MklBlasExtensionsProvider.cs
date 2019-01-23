using IntelMKL.LP64;

//TODO: this should probably call the MKL dll directly, instead of using the package Compute.NET Bindings.
namespace ISAAR.MSolve.LinearAlgebra.Providers.MKL
{
    /// <summary>
    /// Implementation of <see cref="IBlasExtensionsProvider"/> that calls the native dlls of Intel Math Kernel Library. See
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-blas-and-sparse-blas-routines, particularly
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-blas-like-extensions#AEFE7AA3-E2F1-48B1-BCED-92ADB659AA9D
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class MklBlasExtensionsProvider : IBlasExtensionsProvider
    {
        internal static MklBlasExtensionsProvider UniqueInstance { get; } = new MklBlasExtensionsProvider();

        private MklBlasExtensionsProvider() { } // private constructor for singleton pattern

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-axpby#4CF5AEEA-E804-4EF6-BEC0-D2C83CC1DC57
        /// </summary>
        public void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY,
            int incY)
            => Blas.Daxpby(ref n, ref alpha, ref x[offsetX], ref incX, ref beta, ref y[offsetY], ref incY);
    }
}
