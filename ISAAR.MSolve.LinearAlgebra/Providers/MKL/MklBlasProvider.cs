using IntelMKL.LP64;

//TODO: this should probably call the MKL dll directly, instead of using the package Compute.NET Bindings.
namespace ISAAR.MSolve.LinearAlgebra.Providers.MKL
{
    /// <summary>
    /// Implementation of <see cref="IBlasProvider"/> that calls the native dlls of Intel Math Kernel Library. See
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-blas-and-sparse-blas-routines, particularly
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-blas-level-1-routines-and-functions#54BF7621-EE61-46CD-AE5E-A28D786DE838,
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-blas-level-2-routines#708B7246-A41B-4401-80A0-62F4DE57BE79,
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-blas-level-3-routines#1465FCFE-34BB-4C1C-A6E3-7A839CC0842F
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class MklBlasProvider : IBlasProvider
    {
        internal static MklBlasProvider UniqueInstance { get; } = new MklBlasProvider();

        private MklBlasProvider() { } // private constructor for singleton pattern

        #region BLAS Level 1

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-axpy#E25D8E10-0440-4827-BC58-BC71128EA6EE
        /// </summary>
        public void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => Blas.Daxpy(ref n, ref alpha, ref x[offsetX], ref incX, ref y[offsetY], ref incY);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-dot#D4E53C70-D8FA-4095-A800-4203CAFE64FE
        /// </summary>
        public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => Blas.Ddot(ref n, ref x[offsetX], ref incX, ref y[offsetY], ref incY);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-nrm2#EA1DF8E7-FC12-4A82-A804-B62956334C40
        /// </summary>
        public double Dnrm2(int n, double[] x, int offsetX, int incX)
            => Blas.Dnrm2(ref n, ref x[offsetX], ref incX);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-scal#7269DCFE-7235-4690-A69E-D08712D8FC44
        /// </summary>
        public void Dscal(int n, double alpha, double[] x, int offsetX, int incX)
            => Blas.Dscal(ref n, ref alpha, ref x[offsetX], ref incX);
        #endregion

        #region BLAS Level 2

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-gemv#443228C4-626E-48A7-B230-26FB061EACF2
        /// </summary>
        public void Dgemv(TransposeMatrix transA, int m, int n,
            double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
            double beta, double[] y, int offsetY, int incY)
            => Blas.Dgemv(transA.Translate(), ref m, ref n, ref alpha, ref a[offsetA], ref ldA,
                ref x[offsetX], ref incX, ref beta, ref y[offsetY], ref incY);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-spmv#16CB58C4-105B-486C-B6AA-42BB0C721A76
        /// </summary>
        public void Dspmv(StoredTriangle uplo, int n,
            double alpha, double[] a, int offsetA, double[] x, int offsetX, int incX,
            double beta, double[] y, int offsetY, int incY)
            => Blas.Dspmv(uplo.Translate(), ref n, ref alpha, ref a[offsetA],
                ref x[offsetX], ref incX, ref beta, ref y[offsetY], ref incY);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-tpmv#F6666C0E-B843-4E12-9AD4-8898A6EF4018
        /// </summary>
        public void Dtpmv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
            => Blas.Dtpmv(uplo.Translate(), transA.Translate(), diag.Translate(), ref n, ref a[offsetA],
                ref x[offsetX], ref incX);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-tpsv#0EECD264-9871-4097-8AF5-68EEDAE0D00A
        /// </summary>
        public void Dtpsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
            => Blas.Dtpsv(uplo.Translate(), transA.Translate(), diag.Translate(), ref n, ref a[offsetA],
                ref x[offsetX], ref incX);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-trsv#D8733073-F041-4AA1-B82C-123DFA993AD7
        /// </summary>
        public void Dtrsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
            => Blas.Dtrsv(uplo.Translate(), transA.Translate(), diag.Translate(), ref n, ref a[offsetA], ref ldA,
                ref x[offsetX], ref incX);
        #endregion

        #region BLAS Level 3

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-gemm#90EAA001-D4C8-4211-9EA0-B62F5ADE9CF0
        /// </summary>
        public void Dgemm(TransposeMatrix transA, TransposeMatrix transB, int m, int n, int k, double alpha,
            double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC)
            => Blas.Dgemm(transA.Translate(), transB.Translate(), ref m, ref n, ref k, ref alpha, ref a[offsetA], ref ldA,
                ref b[offsetB], ref ldB, ref beta, ref c[offsetC], ref ldC);
        #endregion
    }
}
