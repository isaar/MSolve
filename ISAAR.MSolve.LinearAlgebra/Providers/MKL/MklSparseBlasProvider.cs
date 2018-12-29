using System;
using IntelMKL.LP64;

namespace ISAAR.MSolve.LinearAlgebra.Providers.MKL
{
    /// <summary>
    /// Implementation of <see cref="ISparseBlasProvider"/> that calls the native dlls of Intel Math Kernel Library. See
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-blas-and-sparse-blas-routines, particularly
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-sparse-blas-level-1-routines#08EF100F-3A3A-4CB2-90D3-48B91056BCAD,
    /// https://software.intel.com/en-us/mkl-developer-reference-fortran-sparse-blas-level-2-and-level-3-routines#CCA50EB2-7851-48F7-9CD6-804D7EFD8481
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class MklSparseBlasProvider : ISparseBlasProvider
    {
        internal static MklSparseBlasProvider UniqueInstance { get; } = new MklSparseBlasProvider();

        private MklSparseBlasProvider() { } // private constructor for singleton pattern

        #region Sparse BLAS Level 1

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-axpyi
        /// </summary>
        public void Daxpyi(int nnz, double alpha, double[] valuesX, int[] indicesX, int offsetX, double[] vectorY, int offsetY)
            => CBlas.Daxpyi(nnz, alpha, ref valuesX[offsetX], ref indicesX[offsetX], ref vectorY[offsetY]);

        /// <summary>
        /// See https://software.intel.com/en-us/mkl-developer-reference-fortran-doti
        /// </summary>
        public double Ddoti(int nnz, double[] valuesX, int[] indicesX, int offsetX, double[] vectorY, int offsetY)
            => CBlas.Ddoti(nnz, ref valuesX[offsetX], ref indicesX[offsetX], ref vectorY[offsetY]);
        #endregion

        #region Sparse BLAS Level 2

        public void Dcscgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] colOffsetsA, int[] rowIndicesA,
            double[] vectorX, int offsetX, double[] vectorY, int offsetY)
            => Dcsrgemv(!transposeA, numColsA, numRowsA, valuesA, colOffsetsA, rowIndicesA, vectorX, offsetX, vectorY, offsetY);

        /// <summary>
        /// See
        /// https://software.intel.com/en-us/mkl-developer-reference-fortran-mkl-cspblas-csrgemv#9E1032C5-F844-42D4-A0F0-D62E52D77020
        /// </summary>
        public void Dcsrgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA, int[] colIndicesA,
            double[] vectorX, int offsetX, double[] vectorY, int offsetY)
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
                if (numColsA < numRowsA) // Do not use y.Length, since y can be an unrolled matrix
                {
                    var temp = new double[numRowsA];
                    SpBlas.MklCspblasDcsrgemv("T", ref numRowsA, ref valuesA[0], ref rowOffsetsA[0], ref colIndicesA[0],
                        ref vectorX[offsetX], ref temp[0]);
                    Array.Copy(temp, 0, vectorY, offsetY, numColsA);
                }
                else
                {
                    SpBlas.MklCspblasDcsrgemv("T", ref numRowsA, ref valuesA[0], ref rowOffsetsA[0], ref colIndicesA[0],
                        ref vectorX[offsetX], ref vectorY[offsetY]);
                }
            }
            else
            {
                SpBlas.MklCspblasDcsrgemv("N", ref numRowsA, ref valuesA[0], ref rowOffsetsA[0], ref colIndicesA[0],
                    ref vectorX[offsetX], ref vectorY[offsetY]);
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
        #endregion

        #region Sparse BLAS Level 3

        public void Dcscgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] colOffsetsA,
            int[] rowIndicesA, double[] matrixB, double[] matrixC)
            => Dcsrgemm(!transposeA, numColsA, numColsB, numRowsA, valuesA, colOffsetsA, rowIndicesA, matrixB, matrixC);

        public void Dcsrgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] rowOffsetsA,
            int[] colIndicesA, double[] matrixB, double[] matrixC)
        {
            int ldB = transposeA ? numRowsA : numColsA;
            int ldC = transposeA ? numColsA : numRowsA;

            //TODO: This implementation uses level 2 BLAS. To leverage real lvl 3 BLAS, I would need a different CSR format or  
            //      perhaps the inspector-executor interface.
            for (int j = 0; j < numColsB; ++j)
            {
                Dcsrgemv(transposeA, numRowsA, numColsA, valuesA, rowOffsetsA, colIndicesA, matrixB, j * ldB, matrixC, j * ldC);
            }
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
        #endregion

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
