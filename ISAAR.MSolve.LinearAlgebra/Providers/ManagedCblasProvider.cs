using System;
using DotNumerics.LinearAlgebra.CSLapack;
using static ISAAR.MSolve.LinearAlgebra.Providers.ManagedConstants;

//TODO: find a managed BLAS that supports the methods DotNumerics doesn't or implement them myself (naively).
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public class ManagedCblasProvider : ICblasProvider
    {
        //TODO: perhaps these should not be static.
        private static readonly DAXPY daxpy = new DAXPY();
        private static readonly DDOT ddot = new DDOT();
        private static readonly DGEMM dgemm = new DGEMM();
        private static readonly DGEMV dgemv = new DGEMV();
        private static readonly DNRM2 dnrm2 = new DNRM2();
        private static readonly DSCAL dscal = new DSCAL();

        #region BLAS Level 1
        public void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => daxpy.Run(n, alpha, x, offsetX, incX, ref y, offsetY, incY);

        public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => ddot.Run(n, x, offsetX, incX, y, offsetY, incY);

        public double Dnrm2(int n, double[] x, int offsetX, int incX)
            => dnrm2.Run(n, x, offsetX, incX);

        public void Dscal(int n, double alpha, double[] x, int offsetX, int incX)
            => dscal.Run(n, alpha, ref x, offsetX, incX);
        #endregion

        #region BLAS Level 2

        public void Dgemv(CblasLayout layout, CblasTranspose trans, int m, int n, double alpha, double[] a, int offsetA, int ldA,
            double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
            => dgemv.Run(Translate(layout, trans), m, n, alpha, a, offsetA, ldA, x, offsetX, incX, beta, ref y, offsetY, incY);

        public void Dspmv(CblasLayout layout, CblasTriangular uplo, int n, double alpha, double[] a, int offsetA,
            double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
            => throw new NotImplementedException("The current managed BLAS does not support operations with packed format.");

        public void Dtpmv(CblasLayout layout, CblasTriangular uplo, CblasTranspose trans, CblasDiagonal diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
            => throw new NotImplementedException("The current managed BLAS does not support operations with packed format.");

        public void Dtpsv(CblasLayout layout, CblasTriangular uplo, CblasTranspose trans, CblasDiagonal diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
            => throw new NotImplementedException("The current managed BLAS does not support operations with packed format.");
        #endregion

        #region BLAS Level 3
        public void Dgemm(CblasLayout layout, CblasTranspose transA, CblasTranspose transB, int m, int n, int k, double alpha,
            double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC)
            => dgemm.Run(Translate(layout, transA), Translate(layout, transB), m, n, k, alpha, a, offsetA, ldA, b, offsetB, ldB,
                beta, ref c, offsetC, ldC);
        #endregion

        #region BLAS-like extensions
        public void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, 
            int incY)
        {
            dscal.Run(n, beta, ref y, offsetY, incY);
            daxpy.Run(n, alpha, x, offsetX, incX, ref y, offsetY, incY);
        }
        #endregion
    }
}
