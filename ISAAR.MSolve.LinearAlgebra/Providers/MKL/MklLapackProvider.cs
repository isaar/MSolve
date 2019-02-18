using IntelMKL.LP64;

//TODO: this should probably call the MKL dll directly, instead of using the package Compute.NET Bindings.
namespace ISAAR.MSolve.LinearAlgebra.Providers.MKL
{
    /// <summary>
    /// Implementation of <see cref="ILapackProvider"/> that calls the native dlls of Intel Math Kernel Library. See
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-lapack-routines, particularly
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-lapack-linear-equation-routines,
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-orthogonal-factorizations-lapack-computational-routines#7C288008-E0D0-4E7A-9835-F03B984C3C8E
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class MklLapackProvider : ILapackProvider
    {
        internal static MklLapackProvider UniqueInstance { get; } = new MklLapackProvider();

        private MklLapackProvider() { } // private constructor for singleton pattern

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-gelqf#4E659C89-1230-4052-8454-26AF5066BA46
        /// </summary>
        public void Dgelqf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dgelqf(ref m, ref n, ref a[offsetA], ref ldA, ref tau[offsetTau],
                ref work[offsetWork], ref lWork, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-geqrf#C860486D-506E-44CA-A47B-FA4E5288147D
        /// </summary>
        public void Dgeqrf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, 
            double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dgeqrf(ref m, ref n, ref a[offsetA], ref ldA, ref tau[offsetTau], 
                ref work[offsetWork], ref lWork, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-getrf#42740A2C-4898-4EFA-88B9-94CA6EAAC4DB
        /// </summary>
        public void Dgetrf(int m, int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, ref int info)
            => Lapack.Dgetrf(ref m, ref n, ref a[offsetA], ref ldA, ref ipiv[offsetIpiv], ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-getri#626EB2AE-CA6A-4233-A6FA-04F54EF7A6E6
        /// </summary>
        public void Dgetri(int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv,
            double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dgetri(ref n, ref a[offsetA], ref ldA, ref ipiv[offsetIpiv], ref work[offsetWork], ref lWork, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-getrs#C286AE91-D0E0-44FB-94B4-5C262C037CAF
        /// </summary>
        public void Dgetrs(string transA, int n, int nRhs, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv,
            double[] b, int offsetB, int ldB, ref int info)
            => Lapack.Dgetrs(transA, ref n, ref nRhs, ref a[offsetA], ref ldA, ref ipiv[offsetIpiv], 
                ref b[offsetB], ref ldB, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-orglq#C8D288CD-C41A-46EB-94B6-62C4569FC0DC
        /// </summary>
        public void Dorglq(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, 
            double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dorglq(ref m, ref n, ref k, ref a[offsetA], ref ldA, ref tau[offsetTau],
                ref work[offsetWork], ref lWork, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-orgqr#8CE6CEB0-18AA-4CF6-80F4-D47C9A754019
        /// </summary>
        public void Dorgqr(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau,
            double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dorgqr(ref m, ref n, ref k, ref a[offsetA], ref ldA, ref tau[offsetTau], 
                ref work[offsetWork], ref lWork, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-ormlq#0ECF1EB6-C8B0-4C81-8873-940EFEA2C08B
        /// </summary>
        public void Dormlq(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, 
            int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dormlq(side, transQ, ref m, ref n, ref k, ref a[offsetA], ref ldA, ref tau[offsetTau],
                ref c[offsetC], ref ldC, ref work[offsetWork], ref lWork, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-ormqr#BE4B1C37-0E36-48EB-B44A-A40CACE43821
        /// </summary>
        public void Dormqr(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau,
            int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info)
            => Lapack.Dormqr(side, transQ, ref m, ref n, ref k, ref a[offsetA], ref ldA, ref tau[offsetTau],
                ref c[offsetC], ref ldC, ref work[offsetWork], ref lWork, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-potrf#526C0AD5-B853-4AAC-B27A-E631EE80F066
        /// </summary>
        public void Dpotrf(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
            => Lapack.Dpotrf(uplo, ref n, ref a[offsetA], ref ldA, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-potri#64544865-D81B-4811-84AA-D218A680AA3D
        /// </summary>
        public void Dpotri(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
            => Lapack.Dpotri(uplo, ref n, ref a[offsetA], ref ldA, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-potrs#16D736F0-76B6-4EBA-985F-6EFFD0FB3C41
        /// </summary>
        public void Dpotrs(string uplo, int n, int nRhs, double[] a, int offsetA, int ldA,
            double[] b, int offsetB, int ldB, ref int info)
            => Lapack.Dpotrs(uplo, ref n, ref nRhs, ref a[offsetA], ref ldA, ref b[offsetB], ref ldB, ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-pptrf#A2934477-60D2-40B4-B07D-2AD982989C47
        /// </summary>
        public void Dpptrf(string uplo, int n, double[] a, int offsetA, ref int info)
            => Lapack.Dpptrf(uplo, ref n, ref a[offsetA], ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-pptri#C45F200E-7462-4297-9490-B2433811F035
        /// </summary>
        public void Dpptri(string uplo, int n, double[] a, int offsetA, ref int info)
            => Lapack.Dpptri(uplo, ref n, ref a[offsetA], ref info);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-pptrs#B8D4E277-E28A-42BD-8496-9A72A66EACC3
        /// </summary>
        public void Dpptrs(string uplo, int n, int nRhs, double[] a, int offsetA, double[] b, int offsetB, int ldB, ref int info)
            => Lapack.Dpptrs(uplo, ref n, ref nRhs, ref a[offsetA], ref b[offsetB], ref ldB, ref info);
    }
}
