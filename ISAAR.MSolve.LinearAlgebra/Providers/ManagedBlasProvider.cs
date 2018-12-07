using System;
using DotNumerics.LinearAlgebra.CSLapack;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public class ManagedBlasProvider : IBlasProvider
    {
        #region BLAS Level 1
        //TODO: perhaps these should not be static.
        private static readonly DAXPY daxpy = new DAXPY();
        private static readonly DDOT ddot = new DDOT();
        private static readonly DNRM2 dnrm2 = new DNRM2();
        private static readonly DSCAL dscal = new DSCAL();

        public void Daxpy(int n, double a, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => daxpy.Run(n, a, x, offsetX, incX, ref y, offsetY, incY);

        public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => ddot.Run(n, x, offsetX, incX, y, offsetY, incY);

        public double Dnrm2(int n, double[] x, int offsetX, int incX)
            => dnrm2.Run(n, x, offsetX, incX);

        public void Dscal(int n, double a, double[] x, int offsetX, int incX)
            => dscal.Run(n, a, ref x, offsetX, incX);
        #endregion

        #region BLAS Level 2
        public void FullColMajorTimesVector(int numRows, int numCols, double[] matrix, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }

        public void FullColMajorTransposeTimesVector(int numRows, int numCols, double[] matrix, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }

        public void LowerRowMajorTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }

        public void LowerRowMajorTransposeTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }

        public void SymmColMajorTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }

        public void UpperColMajorTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }

        public void UpperColMajorTransposeTimesVector(int order, double[] matrix, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }
        #endregion

        #region BLAS Level 3
        #endregion

        #region BLAS-like extensions
        public void Daxpby(int n, double a, double[] x, int offsetX, int incX, double b, double[] y, int offsetY, int incY)
        {
            dscal.Run(n, b, ref y, offsetY, incY);
            daxpy.Run(n, a, x, offsetX, incX, ref y, offsetY, incY);
        }
        #endregion
    }
}
