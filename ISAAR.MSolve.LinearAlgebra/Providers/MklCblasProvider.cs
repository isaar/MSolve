using System;
using IntelMKL.LP64;
using static ISAAR.MSolve.LinearAlgebra.Providers.MklConstants;

//TODO: this should probably call the MKL dll directly, instead of using the package Compute.NET Bindings.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public class MklCBlasProvider : ICBlasProvider
    {
        public static MklCBlasProvider UniqueInstance { get; } = new MklCBlasProvider();

        private MklCBlasProvider() { } // private constructor for singleton pattern

        #region BLAS Level 1
        public void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => CBlas.Daxpy(n, alpha, ref x[offsetX], incX, ref y[offsetY], incY);

        public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => CBlas.Ddot(n, ref x[offsetX], incX, ref y[offsetY], incY);

        public double Dnrm2(int n, double[] x, int offsetX, int incX)
            => CBlas.Dnrm2(n, ref x[offsetX], incX);

        public void Dscal(int n, double alpha, double[] x, int offsetX, int incX)
            => CBlas.Dscal(n, alpha, ref x[offsetX], incX);
        #endregion

        #region BLAS Level 2
        public void Dgemv(CBlasLayout layout, CBlasTranspose transA, int m, int n,
            double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX, 
            double beta, double[] y, int offsetY, int incY)
            => CBlas.Dgemv(Translate(layout), Translate(transA), m, n, alpha, ref a[offsetA], ldA, 
                ref x[offsetX], incX, beta, ref y[offsetY], incY);

        public void Dspmv(CBlasLayout layout, CBlasTriangular uplo, int n, 
            double alpha, double[] a, int offsetA, double[] x, int offsetX, int incX,
            double beta, double[] y, int offsetY, int incY)
            => CBlas.Dspmv(Translate(layout), Translate(uplo), n, alpha, ref a[offsetA], 
                ref x[offsetX], incX, beta, ref y[offsetY], incY);

        public void Dtpmv(CBlasLayout layout, CBlasTriangular uplo, CBlasTranspose transA, CBlasDiagonal diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
            => CBlas.Dtpmv(Translate(layout), Translate(uplo), Translate(transA), Translate(diag), n, ref a[offsetA], 
                ref x[offsetX], incX);

        public void Dtpsv(CBlasLayout layout, CBlasTriangular uplo, CBlasTranspose transA, CBlasDiagonal diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
            => CBlas.Dtpsv(Translate(layout), Translate(uplo), Translate(transA), Translate(diag), n, ref a[offsetA],
                ref x[offsetX], incX);

        public void Dtrsv(CBlasLayout layout, CBlasTriangular uplo, CBlasTranspose transA, CBlasDiagonal diag, int n,
            double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
            => CBlas.Dtrsv(Translate(layout), Translate(uplo), Translate(transA), Translate(diag), n, ref a[offsetA], ldA,
                ref x[offsetX], incX);
        #endregion

        #region BLAS Level 3
        public void Dgemm(CBlasLayout layout, CBlasTranspose transA, CBlasTranspose transB, int m, int n, int k, double alpha,
            double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC)
            => CBlas.Dgemm(Translate(layout), Translate(transA), Translate(transB), m, n, k, alpha, ref a[offsetA], ldA,
                ref b[offsetB], ldB, beta, ref c[offsetC], ldC);
        #endregion

        #region BLAS-like extensions
        public void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, 
            int incY)
            => CBlas.Daxpby(n, alpha, ref x[offsetX], incX, beta, ref y[offsetY], incY);
        #endregion

        
    }
}
