using IntelMKL.LP64;

//TODO: this should probably call the MKL dll directly, instead of using the package Compute.NET Bindings.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    internal class MklLapackProvider : ILapackProvider
    {
        internal static MklLapackProvider UniqueInstance { get; } = new MklLapackProvider();

        private MklLapackProvider() { } // private constructor for singleton pattern

        public void Dgelqf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dgelqf(ref m, ref n, ref a[offsetA], ref ldA, ref tau[offsetTau],
                ref work[offsetWork], ref lWork, ref info);

        public void Dgeqrf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, 
            double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dgeqrf(ref m, ref n, ref a[offsetA], ref ldA, ref tau[offsetTau], 
                ref work[offsetWork], ref lWork, ref info);

        public void Dgetrf(int m, int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, ref int info)
            => Lapack.Dgetrf(ref m, ref n, ref a[offsetA], ref ldA, ref ipiv[offsetIpiv], ref info);

        public void Dgetri(int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv,
            double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dgetri(ref n, ref a[offsetA], ref ldA, ref ipiv[offsetIpiv], ref work[offsetWork], ref lWork, ref info);

        public void Dgetrs(string transA, int n, int nRhs, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv,
            double[] b, int offsetB, int ldB, ref int info)
            => Lapack.Dgetrs(transA, ref n, ref nRhs, ref a[offsetA], ref ldA, ref ipiv[offsetIpiv], 
                ref b[offsetB], ref ldB, ref info);

        public void Dorglq(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, 
            double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dorglq(ref m, ref n, ref k, ref a[offsetA], ref ldA, ref tau[offsetTau],
                ref work[offsetWork], ref lWork, ref info);

        public void Dorgqr(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau,
            double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dorgqr(ref m, ref n, ref k, ref a[offsetA], ref ldA, ref tau[offsetTau], 
                ref work[offsetWork], ref lWork, ref info);

        public void Dormlq(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, 
            int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dormlq(side, transQ, ref m, ref n, ref k, ref a[offsetA], ref ldA, ref tau[offsetTau],
                ref c[offsetC], ref ldC, ref work[offsetWork], ref lWork, ref info);

        public void Dormqr(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau,
            int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dormqr(side, transQ, ref m, ref n, ref k, ref a[offsetA], ref ldA, ref tau[offsetTau],
                ref c[offsetC], ref ldC, ref work[offsetWork], ref lWork, ref info);

        public void Dpotrf(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
            => Lapack.Dpotrf(uplo, ref n, ref a[offsetA], ref ldA, ref info);

        public void Dpotri(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
            => Lapack.Dpotri(uplo, ref n, ref a[offsetA], ref ldA, ref info);

        public void Dpotrs(string uplo, int n, int nRhs, double[] a, int offsetA, int ldA,
            double[] b, int offsetB, int ldB, ref int info)
            => Lapack.Dpotrs(uplo, ref n, ref nRhs, ref a[offsetA], ref ldA, ref b[offsetB], ref ldB, ref info);

        public void Dpptrf(string uplo, int n, double[] a, int offsetA, ref int info)
            => Lapack.Dpptrf(uplo, ref n, ref a[offsetA], ref info);

        public void Dpptri(string uplo, int n, double[] a, int offsetA, ref int info)
            => Lapack.Dpptri(uplo, ref n, ref a[offsetA], ref info);

        public void Dpptrs(string uplo, int n, int nRhs, double[] a, int offsetA, double[] b, int offsetB, int ldB, ref int info)
            => Lapack.Dpptrs(uplo, ref n, ref nRhs, ref a[offsetA], ref b[offsetB], ref ldB, ref info);
    }
}
