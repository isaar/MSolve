using IntelMKL.LP64;

//TODO: this should probably call the MKL dll directly, instead of using the package Compute.NET Bindings.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public class MklBlasProvider : IBlasProvider
    {
        public static MklBlasProvider UniqueInstance { get; } = new MklBlasProvider();

        private MklBlasProvider() { } // private constructor for singleton pattern

        #region BLAS Level 1
        public void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => Blas.Daxpy(ref n, ref alpha, ref x[offsetX], ref incX, ref y[offsetY], ref incY);

        public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => Blas.Ddot(ref n, ref x[offsetX], ref incX, ref y[offsetY], ref incY);

        public double Dnrm2(int n, double[] x, int offsetX, int incX)
            => Blas.Dnrm2(ref n, ref x[offsetX], ref incX);

        public void Dscal(int n, double alpha, double[] x, int offsetX, int incX)
            => Blas.Dscal(ref n, ref alpha, ref x[offsetX], ref incX);
        #endregion

        #region BLAS Level 2
        public void Dgemv(TransposeMatrix transA, int m, int n,
            double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
            double beta, double[] y, int offsetY, int incY)
            => Blas.Dgemv(transA.Translate(), ref m, ref n, ref alpha, ref a[offsetA], ref ldA,
                ref x[offsetX], ref incX, ref beta, ref y[offsetY], ref incY);

        public void Dspmv(StoredTriangle uplo, int n,
            double alpha, double[] a, int offsetA, double[] x, int offsetX, int incX,
            double beta, double[] y, int offsetY, int incY)
            => Blas.Dspmv(uplo.Translate(), ref n, ref alpha, ref a[offsetA],
                ref x[offsetX], ref incX, ref beta, ref y[offsetY], ref incY);

        public void Dtpmv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
            => Blas.Dtpmv(uplo.Translate(), transA.Translate(), diag.Translate(), ref n, ref a[offsetA],
                ref x[offsetX], ref incX);

        public void Dtpsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
            => Blas.Dtpsv(uplo.Translate(), transA.Translate(), diag.Translate(), ref n, ref a[offsetA],
                ref x[offsetX], ref incX);

        public void Dtrsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
            => Blas.Dtrsv(uplo.Translate(), transA.Translate(), diag.Translate(), ref n, ref a[offsetA], ref ldA,
                ref x[offsetX], ref incX);
        #endregion

        #region BLAS Level 3
        public void Dgemm(TransposeMatrix transA, TransposeMatrix transB, int m, int n, int k, double alpha,
            double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC)
            => Blas.Dgemm(transA.Translate(), transB.Translate(), ref m, ref n, ref k, ref alpha, ref a[offsetA], ref ldA,
                ref b[offsetB], ref ldB, ref beta, ref c[offsetC], ref ldC);
        #endregion
    }
}
