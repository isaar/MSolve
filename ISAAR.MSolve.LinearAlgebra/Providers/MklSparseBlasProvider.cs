using System;
using System.Collections.Generic;
using System.Text;
using IntelMKL.LP64;

//TODO: Should I use transpose flags instead of providing different methods. The caller would have to repeat checking this
//      flag (granted, it is not that costly) and I would need to write a lot of nested if...else... clauses.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Delegates BLAS operations to the highly optimized native dlls provided by Intel MKL.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class MklSparseBlasProvider : ISparseBlasProvider
    {
        public void CscTimesVector(int numCols, double[] values, int[] colOffsets, int[] rowIndices, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }

        public void CscTransposeTimesVector(int numCols, double[] values, int[] colOffsets, int[] rowIndices, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }

        public void CsrTimesVector(int numRows, double[] values, int[] rowOffsets, int[] colIndices, double[] x, double[] y)
        {
            SpBlas.MklCspblasDcsrgemv("N", ref numRows, ref values[0], ref rowOffsets[0], ref colIndices[0],
                    ref x[0], ref y[0]);
        }

        public void CsrTransposeTimesVector(int numRows, double[] values, int[] rowOffsets, int[] colIndices, double[] x, 
            double[] y)
        {
            // SparseBLAS overwrites memory equal to the number of matrix rows. If the rhs vector is shorter than that, the 
            // remaining entries are overwritten with 0. However, this messes up the managed vector objects, since 
            // important data is overwritten (why doesn't it throw access violation exception?). For now I am going to 
            // use a temp array and then copy the relevant part. This problem does not seem to appear in the untransposed 
            // version of the method.
            //TODO: Try using the SparseBLAS inspector-executor routines, instead of the deprecated dcsrgemv().
            if (y.Length < numRows)
            {
                var temp = new double[numRows];
                SpBlas.MklCspblasDcsrgemv("T", ref numRows, ref values[0], ref rowOffsets[0], ref colIndices[0],
                    ref x[0], ref temp[0]);
                Array.Copy(temp, y, y.Length);
            }
            else
            {
                SpBlas.MklCspblasDcsrgemv("T", ref numRows, ref values[0], ref rowOffsets[0], ref colIndices[0],
                    ref x[0], ref y[0]);
            }
        }

        public void SkylineTimesVector(int order, double[] values, int[] diagOffsets, double[] x, double[] y)
        {
            throw new NotImplementedException();
        }
    }
}
