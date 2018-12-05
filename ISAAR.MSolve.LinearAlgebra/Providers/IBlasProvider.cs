using System;
using System.Collections.Generic;
using System.Text;

//TODO: this and its implementations should be internal. The user should select a Provider that will then specify different BLAS,
//      SparseBLAS, LAPACK, etc. providers. The implementations shoudl also be singletons or enums.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public interface IBlasProvider
    {
        void FullColMajorTimesVector(int numRows, int numCols, double[] matrix, double[] x, double[] y);
        void FullColMajorTransposeTimesVector(int numRows, int numCols, double[] matrix, double[] x, double[] y);
        void LowerRowMajorTimesVector(int order, double[] matrix, double[] x, double[] y);
        void LowerRowMajorTransposeTimesVector(int order, double[] matrix, double[] x, double[] y);
        void SymmColMajorTimesVector(int order, double[] matrix, double[] x, double[] y);
        void UpperColMajorTimesVector(int order, double[] matrix, double[] x, double[] y);
        void UpperColMajorTransposeTimesVector(int order, double[] matrix, double[] x, double[] y);
    }
}
