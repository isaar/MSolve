using DotNumerics.LinearAlgebra.CSLapack;

namespace ISAAR.MSolve.LinearAlgebra.Providers.Managed
{
    /// <summary>
    /// Provides managed C# implementations of the linear algebra operations defined by <see cref="IBlasExtensionsProvider"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class ManagedBlasExtensionsProvider : IBlasExtensionsProvider
    {
        //TODO: perhaps these should not be static.
        private static readonly DAXPY daxpy = new DAXPY();
        private static readonly DSCAL dscal = new DSCAL();

        internal static ManagedBlasExtensionsProvider UniqueInstance { get; } = new ManagedBlasExtensionsProvider();

        private ManagedBlasExtensionsProvider() { } // private constructor for singleton pattern

        public void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
        {
            dscal.Run(n, beta, ref y, offsetY, incY);
            daxpy.Run(n, alpha, x, offsetX, incX, ref y, offsetY, incY);
        }
    }
}
