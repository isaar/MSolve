using System;
using DotNumerics.LinearAlgebra.CSLapack;
using ISAAR.MSolve.LinearAlgebra.Providers.Implementations;
using static ISAAR.MSolve.LinearAlgebra.Providers.ManagedConstants;

//TODO: find a managed BLAS that supports the methods DotNumerics doesn't.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public class ManagedCBlasProvider : ICBlasProvider
    {
        //TODO: perhaps these should not be static.
        private static readonly DAXPY daxpy = new DAXPY();
        private static readonly DDOT ddot = new DDOT();
        private static readonly DGEMM dgemm = new DGEMM();
        private static readonly DGEMV dgemv = new DGEMV();
        private static readonly DNRM2 dnrm2 = new DNRM2();
        private static readonly DSCAL dscal = new DSCAL();
        private static readonly DTRSV dtrsv = new DTRSV();

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

        public void Dgemv(CBlasLayout layout, CBlasTranspose transA, int m, int n,
            double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX, 
            double beta, double[] y, int offsetY, int incY)
            => dgemv.Run(Translate(layout, transA), m, n, alpha, a, offsetA, ldA, x, offsetX, incX, beta, ref y, offsetY, incY);

        public void Dspmv(CBlasLayout layout, CBlasTriangular uplo, int n, 
            double alpha, double[] a, int offsetA, double[] x, int offsetX, int incX, 
            double beta, double[] y, int offsetY, int incY)
        {
            // y = alpha * L * x + beta * y 
            CblasLevel2Implementations.LowerTimesVectorPackedRowMajor(
                CblasLevel2Implementations.Diagonal.Regular, n, alpha, a, offsetA, x, offsetX, incX, beta, y, offsetY, incY);

            // y = alpha * U * x + y, where U has 0 diagonal
            CblasLevel2Implementations.UpperTimesVectorPackedColMajor(
                CblasLevel2Implementations.Diagonal.Zero, n, alpha, a, offsetA, x, offsetX, incX, 1.0, y, offsetY, incY);
        }

        public void Dtpmv(CBlasLayout layout, CBlasTriangular uplo, CBlasTranspose transA, CBlasDiagonal diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
        {
            // The copy may be avoidable in trangular operations, if we start the dot products from the bottom
            var input = new double[x.Length];
            Array.Copy(x, input, x.Length);

            CblasLevel2Implementations.Diagonal managedDiag = (diag == CBlasDiagonal.NonUnit) ? 
                CblasLevel2Implementations.Diagonal.Regular : CblasLevel2Implementations.Diagonal.Unit;
            if (UseUpperImplementation(uplo, transA, layout))
            {
                CblasLevel2Implementations.UpperTimesVectorPackedColMajor(
                    managedDiag, n, 1.0, a, offsetA, input, offsetX, incX, 0.0, x, offsetX, incX);
            }
            else
            {
                CblasLevel2Implementations.LowerTimesVectorPackedRowMajor(
                    managedDiag, n, 1.0, a, offsetA, input, offsetX, incX, 0.0, x, offsetX, incX);
            }
        }

        public void Dtpsv(CBlasLayout layout, CBlasTriangular uplo, CBlasTranspose transA, CBlasDiagonal diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
        {
            bool unit = (diag == CBlasDiagonal.Unit) ? true : false;
            if (UseUpperImplementation(uplo, transA, layout))
            {
                CblasLevel2Implementations.BackSubstitutionPackedColMajor(unit, n, a, offsetA, x, offsetX, incX);
            }
            else CblasLevel2Implementations.ForwardSubstitutionPackedRowMajor(unit, n, a, offsetA, x, offsetX, incX);
        }

        public void Dtrsv(CBlasLayout layout, CBlasTriangular uplo, CBlasTranspose transA, CBlasDiagonal diag, int n,
            double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
            => dtrsv.Run(Translate(uplo), Translate(layout, transA), Translate(diag), n, a, offsetA, ldA, ref x, offsetX, incX);
        #endregion

        #region BLAS Level 3
        public void Dgemm(CBlasLayout layout, CBlasTranspose transA, CBlasTranspose transB, int m, int n, int k, double alpha,
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

        //TODO: I have not tested all combinations
        private static bool UseUpperImplementation(CBlasTriangular uplo, CBlasTranspose transA, CBlasLayout layout)
        {
            if (transA == CBlasTranspose.ConjugateTranspose)
                throw new ArgumentException("Cannot use conjugate transpose operations for double matrices and vectors.");

            // Upper triangular combinations
            if (uplo == CBlasTriangular.Upper && transA == CBlasTranspose.NoTranspose && layout == CBlasLayout.ColMajor)
            {
                return true;
            }
            if (uplo == CBlasTriangular.Upper && transA == CBlasTranspose.NoTranspose && layout == CBlasLayout.RowMajor)
            {
                return false;
            }
            if (uplo == CBlasTriangular.Upper && transA == CBlasTranspose.Transpose && layout == CBlasLayout.ColMajor)
            {
                return false;
            }
            if (uplo == CBlasTriangular.Upper && transA == CBlasTranspose.Transpose && layout == CBlasLayout.RowMajor)
            {
                return true;
            }

            // Lower triangular combinations
            if (uplo == CBlasTriangular.Lower && transA == CBlasTranspose.NoTranspose && layout == CBlasLayout.RowMajor)
            {
                return false;
            }
            if (uplo == CBlasTriangular.Lower && transA == CBlasTranspose.NoTranspose && layout == CBlasLayout.ColMajor)
            {
                return true;
            }
            if (uplo == CBlasTriangular.Lower && transA == CBlasTranspose.Transpose && layout == CBlasLayout.RowMajor)
            {
                return true;
            }
            if (uplo == CBlasTriangular.Lower && transA == CBlasTranspose.Transpose && layout == CBlasLayout.ColMajor)
            {
                return false;
            }

            throw new Exception("This code should not have been reached");
        }
    }
}
