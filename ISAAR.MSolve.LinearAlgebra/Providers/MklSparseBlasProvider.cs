using System;
using System.Collections.Generic;
using System.Text;
using IntelMKL.LP64;

//TODO: Should I use transposeA flags instead of providing different methods. The caller would have to repeat checking this
//      flag (granted, it is not that costly) and I would need to write a lot of nested if...else... clauses.
//TODO: This should be a singleton.
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

        public void Dcscgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] colOffsetsA,
            int[] rowIndicesA, double[] b, double[] c)
            => Dcsrgemm(!transposeA, numColsA, numColsB, numRowsA, valuesA, colOffsetsA, rowIndicesA, b, c);

        public void Dcscgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] colOffsetsA, int[] rowIndicesA,
            double[] x, int offsetX, double[] y, int offsetY)
            => Dcsrgemv(!transposeA, numColsA, numRowsA, valuesA, colOffsetsA, rowIndicesA, x, offsetX, y, offsetY);
        
        public void Dcsrgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] rowOffsetsA, 
            int[] colIndicesA, double[] b, double[] c)
        {
            int ldB = transposeA ? numRowsA : numColsA;
            int ldC = transposeA ? numColsA : numRowsA;

            //TODO: This implementation uses level 2 BLAS. To leverage real lvl 3 BLAS, I would need a different CSR format or  
            //      perhaps the inspector-executor interface.
            for (int j = 0; j < numColsB; ++j)
            {
                Dcsrgemv(transposeA, numRowsA, numColsA, valuesA, rowOffsetsA, colIndicesA, b, j * ldB, c, j * ldC);
            }
        }

        /// <summary>
        /// See
        /// https://software.intel.com/en-us/mkl-developer-reference-c-mkl-cspblas-csrgemv#D840F0E5-E41A-4E91-94D2-FEB320F93E91
        /// </summary>
        public void Dcsrgemv(bool transposeA, int numRowsA, int numCols, double[] valuesA, int[] rowOffsetsA, int[] colIndicesA, 
            double[] x, int offsetX, double[] y, int offsetY)
        {
            if (transposeA)
            {
                // This MKL function is only for square CSR matrices. We can use it for rectangular too, but it overwrites memory 
                // of the rhs vector y equal to the number of matrix rows. If the rhs vector is shorter than that, then the 
                // remaining entries are overwritten with 0. However, this messes up the managed vector objects, since 
                // important data is overwritten (why doesn't it throw access violation exception?). For now I am going to 
                // use a temp array and then copy the relevant part. This problem does not seem to appear in the untransposeAd 
                // version of the method.
                //TODO: Try using the SparseBLAS inspector-executor routines, instead of the deprecated dcsrgemv().
                if (numCols < numRowsA) // Do not use y.Length, since y can be an unrolled matrix
                {
                    var temp = new double[numRowsA];
                    SpBlas.MklCspblasDcsrgemv("T", ref numRowsA, ref valuesA[0], ref rowOffsetsA[0], ref colIndicesA[0],
                        ref x[offsetX], ref temp[0]);
                    Array.Copy(temp, 0, y, offsetY, numCols);
                }
                else
                {
                    SpBlas.MklCspblasDcsrgemv("T", ref numRowsA, ref valuesA[0], ref rowOffsetsA[0], ref colIndicesA[0],
                        ref x[offsetX], ref y[offsetY]);
                }
            }
            else
            {
                SpBlas.MklCspblasDcsrgemv("N", ref numRowsA, ref valuesA[0], ref rowOffsetsA[0], ref colIndicesA[0],
                    ref x[offsetX], ref y[offsetY]);
            }
        }

        /// <summary>
        /// Matrix vector multiplication y = A * x, with A being a triangular matrix in Skyline format. See
        /// https://software.intel.com/en-us/mkl-developer-reference-c-mkl-skymv#7C27501E-3893-4443-8257-0E5F4E905C29
        /// Warning: Intel MKL's Skyline format uses 1-based indexing and the non zero entries of each column are ordered from
        /// the top to the diagonal, if the upper triangle is stored. See
        /// https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-skyline-matrix-storage-format
        /// In constrast, in the current version of our Skyline matrix, 0-based indexing is used and the entries of each column 
        /// are ordered from the diagonal to the top.
        /// </summary>
        public void Dskymv(bool upper, int orderA, double[] valuesA, int[] diagOffsetsA, double[] x, double[] y)
        {
            string trans = "N"; //TODO: not sure about this
            string matdscrA = "SUNF"; // symmetric, upper, non-unit, fortran indexing
            double alpha = 1.0;
            double beta = 0.0;
            SpBlas.MklDskymv(trans, ref orderA, ref orderA, ref alpha, matdscrA, ref valuesA[0], ref diagOffsetsA[0],
                ref x[0], ref beta, ref y[0]);
        }

        /// <summary>
        /// Solve multiple linear systems X = inv(A) * B, with A being a triangular matrix in Skyline format. See
        /// https://software.intel.com/node/0a9d506f-d424-4651-8e68-16625ed412e7#0A9D506F-D424-4651-8E68-16625ED412E7
        /// Warning: Intel MKL's Skyline format uses 1-based indexing and the non zero entries of each column are ordered from
        /// the top to the diagonal, if the upper triangle is stored. See
        /// https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-skyline-matrix-storage-format
        /// In constrast, in the current version of our Skyline matrix, 0-based indexing is used and the entries of each column 
        /// are ordered from the diagonal to the top.
        /// </summary>
        public void Dskysm(bool upper, int n, int nRhs, double[] valuesA, int[] colOffsetsA, double[] b, int ldB,
            double[] x, int ldX)
        {
            string trans = "N"; //TODO: not sure about this
            string matdscrA = "SUNF"; // symmetric, upper, non-unit, fortran indexing
            double alpha = 1.0;
            SpBlas.MklDskysm(trans, ref n, ref nRhs, ref alpha, matdscrA, ref valuesA[0], ref colOffsetsA[0],
                ref b[0], ref ldB, ref x[0], ref ldX);
        }

        /// <summary>
        /// Linear system solution x = inv(A) * b, with A being a triangular matrix in Skyline format. See
        /// https://software.intel.com/node/0a9d506f-d424-4651-8e68-16625ed412e7#0A9D506F-D424-4651-8E68-16625ED412E7
        /// Warning: Intel MKL's Skyline format uses 1-based indexing and the non zero entries of each column are ordered from
        /// the top to the diagonal, if the upper triangle is stored. See
        /// https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-skyline-matrix-storage-format
        /// In constrast, in the current version of our Skyline matrix, 0-based indexing is used and the entries of each column 
        /// are ordered from the diagonal to the top.
        /// </summary>
        public void Dskysv(bool upper, int orderA, double[] valuesA, int[] colOffsetsA, double[] b, double[] x)
        {
            string trans = "N"; //TODO: not sure about this
            string matdscrA = "SUNF"; // symmetric, upper, non-unit, fortran indexing
            double alpha = 1.0;
            SpBlas.MklDskysv(trans, ref orderA, ref alpha, matdscrA, ref valuesA[0], ref colOffsetsA[0],
                ref b[0], ref x[0]);
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
