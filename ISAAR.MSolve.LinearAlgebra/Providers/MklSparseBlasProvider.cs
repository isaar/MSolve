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
        public void Daxpyi(int nnz, double[] alpha, double[] x, int offsetX, int[] indicesX, double[] y, int offsetY)
        {
            throw new NotImplementedException();
        }

        public void Dcscgemv(bool transpose, int numRows, int numCols, double[] values, int[] colOffsets, int[] rowIndices,
            double[] x, double[] y)
            => Dcsrgemv(!transpose, numCols, numRows, values, colOffsets, rowIndices, x, y);

        /// <summary>
        /// See
        /// https://software.intel.com/en-us/mkl-developer-reference-c-mkl-cspblas-csrgemv#D840F0E5-E41A-4E91-94D2-FEB320F93E91
        /// </summary>
        public void Dcsrgemv(bool transpose, int numRows, int numColumns,double[] values, int[] rowOffsets, int[] colIndices, 
            double[] x, double[] y)
        {
            if (transpose)
            {
                // This MKL function is only for square CSR matrices. We can use it for rectangular too, but it overwrites memory 
                // of the rhs vector y equal to the number of matrix rows. If the rhs vector is shorter than that, then the 
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
            else
            {
                SpBlas.MklCspblasDcsrgemv("N", ref numRows, ref values[0], ref rowOffsets[0], ref colIndices[0],
                    ref x[0], ref y[0]);
            }
        }

        /// <summary>
        /// Matrix vector multiplication, with a triangular matrix in Skyline format. See
        /// https://software.intel.com/en-us/mkl-developer-reference-c-mkl-skymv#7C27501E-3893-4443-8257-0E5F4E905C29
        /// Warning: Intel MKL's Skyline format uses 1-based indexing and the non zero entries of each column are ordered from
        /// the top to the diagonal, if the upper triangle is stored. See
        /// https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-skyline-matrix-storage-format
        /// In constrast, in the current version of our Skyline matrix, 0-based indexing is used and the entries of each column 
        /// are ordered from the diagonal to the top.
        /// </summary>
        public void Dskymv(bool upper, int order, double[] values, int[] diagOffsets, double[] lhs, double[] rhs)
        {
            //int[] diagOffsets1Based = Convert0To1BasedIndexing(diagOffsets); // This method accepts 1-based indexing.
            //SpBlas.MklDskymv();
            throw new NotImplementedException();
        }

        /// <summary>
        /// Linear system solution, with a triangular matrix in Skyline format. See
        /// https://software.intel.com/node/0a9d506f-d424-4651-8e68-16625ed412e7#0A9D506F-D424-4651-8E68-16625ED412E7
        /// Warning: Intel MKL's Skyline format uses 1-based indexing and the non zero entries of each column are ordered from
        /// the top to the diagonal, if the upper triangle is stored. See
        /// https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-skyline-matrix-storage-format
        /// In constrast, in the current version of our Skyline matrix, 0-based indexing is used and the entries of each column 
        /// are ordered from the diagonal to the top.
        /// </summary>
        public void Dskysv(bool upper, int order, double[] values, int[] diagOffsets, double[] x, double[] y)
        {
            //int[] diagOffsets1Based = Convert0To1BasedIndexing(diagOffsets); // This method accepts 1-based indexing.
            //SpBlas.MklDskysv();
            throw new NotImplementedException();
        }

        //TODO: perhaps the 1-based indexing array can be stored in the matrix or cached here instead of rebuilding it 
        //      each time.
        private static int[] Convert0To1BasedIndexing(int[] zeroBasedIndices)
        {
            int count = zeroBasedIndices.Length;
            var oneBasesIndices = new int[count];
            for (int i = 0; i < count; ++i) oneBasesIndices[i] = zeroBasedIndices[i] + 1;
            return oneBasesIndices;
        }
    }
}
