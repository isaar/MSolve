//TODO: this and its implementations should be internal. The user should select a Provider that will then specify different BLAS,
//      SparseBLAS, LAPACK, etc. providers. The implementations shoudl also be singletons or enums.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public interface IBlasProvider
    {
        #region BLAS Level 1
        void Daxpy(int n, double a, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY);
        double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY);
        double Dnrm2(int n, double[] x, int offsetX, int incX);
        void Dscal(int n, double a, double[] x, int offsetX, int incX);
        #endregion

        #region BLAS Level 2
        //TODO: Use the regular BLAS names and parameters for these.
        void FullColMajorTimesVector(int numRows, int numCols, double[] matrix, double[] x, double[] y);
        void FullColMajorTransposeTimesVector(int numRows, int numCols, double[] matrix, double[] x, double[] y);
        void LowerRowMajorTimesVector(int order, double[] matrix, double[] x, double[] y);
        void LowerRowMajorTransposeTimesVector(int order, double[] matrix, double[] x, double[] y);
        void SymmColMajorTimesVector(int order, double[] matrix, double[] x, double[] y);
        void UpperColMajorTimesVector(int order, double[] matrix, double[] x, double[] y);
        void UpperColMajorTransposeTimesVector(int order, double[] matrix, double[] x, double[] y);
        #endregion

        #region BLAS Level 3
        #endregion

        #region BLAS-like extensions
        void Daxpby(int n, double a, double[] x, int offsetX, int incX, double b, double[] y, int offsetY, int incY);
        #endregion
    }
}
