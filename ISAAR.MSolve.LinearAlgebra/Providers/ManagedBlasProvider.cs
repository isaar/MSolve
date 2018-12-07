using System;
using System.Collections.Generic;
using System.Text;
using DotNumerics.LinearAlgebra.CSLapack;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public class ManagedBlasProvider : IBlasProvider
    {
        private static readonly DAXPY daxpy = new DAXPY();

        public void Daxpy(int n, double a, double[] x, int offsetX, int incX, double[] y, int offsetY, int incy)
            => daxpy.Run(n, a, x, offsetX, incX, ref y, offsetY, incy);

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
    }
}
