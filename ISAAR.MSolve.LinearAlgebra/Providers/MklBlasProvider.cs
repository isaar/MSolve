using System;
using System.Collections.Generic;
using System.Text;
using IntelMKL.LP64;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public class MklBlasProvider : IBlasProvider
    {
        public void Daxpy(int n, double a, double[] x, int offsetX, int incX, double[] y, int offsetY, int incy)
            => CBlas.Daxpy(n, a, ref x[offsetX], incX, ref y[offsetY], 1);

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
    }
}
