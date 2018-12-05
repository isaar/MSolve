using System;
using System.Collections.Generic;
using System.Text;

//TODO: this and its implementations should be internal. The user should select a Provider that will then specify different BLAS,
//      SparseBLAS, LAPACK, etc. providers. The implementations shoudl also be singletons or enums.
//TODO: some dimensions are redundant, since they can be read from the indexing arrays.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public interface ISparseBlasProvider
    {
        void CscTimesVector(int numCols, double[] values, int[] colOffsets, int[] rowIndices, double[] x, double[] y);
        void CscTransposeTimesVector(int numCols, double[] values, int[] colOffsets, int[] rowIndices, double[] x, double[] y);
        void CsrTimesVector(int numRows, double[] values, int[] rowOffsets, int[] colIndices, double[] x, double[] y);
        void CsrTransposeTimesVector(int numRows, double[] values, int[] rowOffsets, int[] colIndices, double[] x, double[] y);
        void SkylineTimesVector(int order, double[] values, int[] diagOffsets, double[] x, double[] y);
    }
}
