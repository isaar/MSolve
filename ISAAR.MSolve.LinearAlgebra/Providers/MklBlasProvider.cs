using System;
using IntelMKL.LP64;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public class MklBlasProvider : IBlasProvider
    {
        #region BLAS Level 1
        public void Daxpy(int n, double a, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => CBlas.Daxpy(n, a, ref x[offsetX], incX, ref y[offsetY], incY);

        public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => CBlas.Ddot(n, ref x[offsetX], incX, ref y[offsetY], incY);

        public double Dnrm2(int n, double[] x, int offsetX, int incX)
            => CBlas.Dnrm2(n, ref x[offsetX], incX);

        public void Dscal(int n, double a, double[] x, int offsetX, int incX)
            => CBlas.Dscal(n, a, ref x[offsetX], incX);
        #endregion

        #region BLAS Level 2
        public void FullColMajorTimesVector(int numRows, int numCols, double[] matrix, double[] x, double[] y)
        {
            CBlas.Dgemv(CBLAS_LAYOUT.CblasColMajor, CBLAS_TRANSPOSE.CblasNoTrans, numRows, numCols,
                1.0, ref matrix[0], numRows,
                ref x[0], 1,
                0.0, ref y[0], 1);
        }

        public void FullColMajorTransposeTimesVector(int numRows, int numCols, double[] matrix, double[] x, double[] y)
        {
            CBlas.Dgemv(CBLAS_LAYOUT.CblasColMajor, CBLAS_TRANSPOSE.CblasTrans, numRows, numCols,
                1.0, ref matrix[0], numRows,
                ref x[0], 1,
                0.0, ref y[0], 1);
        }

        public void LowerRowMajorTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            Array.Copy(x, y, order);
            CBlas.Dtpmv(CBLAS_LAYOUT.CblasRowMajor, CBLAS_UPLO.CblasLower, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit, 
                order, ref matrix[0], ref y[0], 1);
        }

        public void LowerRowMajorTransposeTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            Array.Copy(x, y, order);
            CBlas.Dtpmv(CBLAS_LAYOUT.CblasRowMajor, CBLAS_UPLO.CblasLower, CBLAS_TRANSPOSE.CblasTrans, CBLAS_DIAG.CblasNonUnit,
                order, ref matrix[0], ref y[0], 1);
        }

        public void SymmColMajorTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            CBlas.Dspmv(CBLAS_LAYOUT.CblasColMajor, CBLAS_UPLO.CblasUpper, order,
                1.0, ref matrix[0], ref x[0], 1, 0.0, ref y[0], 1);
        }

        public void UpperColMajorTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            Array.Copy(x, y, order);
            CBlas.Dtpmv(CBLAS_LAYOUT.CblasColMajor, CBLAS_UPLO.CblasUpper, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit,
                order, ref matrix[0], ref y[0], 1);
        }

        public void UpperColMajorTransposeTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            Array.Copy(x, y, order);
            CBlas.Dtpmv(CBLAS_LAYOUT.CblasColMajor, CBLAS_UPLO.CblasUpper, CBLAS_TRANSPOSE.CblasTrans, CBLAS_DIAG.CblasNonUnit,
                order, ref matrix[0], ref y[0], 1);
        }
        #endregion

        #region BLAS Level 3
        #endregion

        #region BLAS-like extensions
        public void Daxpby(int n, double a, double[] x, int offsetX, int incX, double b, double[] y, int offsetY, int incY)
            => CBlas.Daxpby(n, a, ref x[offsetX], incX, b, ref y[offsetY], incY);
        #endregion
    }
}
