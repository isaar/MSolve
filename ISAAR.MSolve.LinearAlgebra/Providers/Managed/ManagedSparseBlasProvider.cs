//TODO: At some point I should change my Skyline format to match the one used in MKL. Then the SparseBLAS operations can be 
//      interchangeable.
namespace ISAAR.MSolve.LinearAlgebra.Providers.Managed
{
    /// <summary>
    /// Provides managed C# implementations of the linear algebra operations defined by <see cref="IBlasProvider"/>.
    /// Uses custom C# code (usually unoptimized).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class ManagedSparseBlasProvider : ISparseBlasProvider
    {
        internal static ManagedSparseBlasProvider UniqueInstance { get; } = new ManagedSparseBlasProvider();

        private ManagedSparseBlasProvider() { } // private constructor for singleton pattern

        #region Sparse BLAS Level 1

        public void Daxpyi(int nnz, double alpha, double[] valuesX, int[] indicesX, int offsetX, double[] vectorY, int offsetY)
            => SparseBlasImplementations.AlphaTimesSparsePlusDenseVector(nnz, alpha, valuesX, indicesX, offsetX,
                vectorY, offsetY);

        public double Ddoti(int nnz, double[] valuesX, int[] indicesX, int offsetX, double[] vectorY, int offsetY)
            => SparseBlasImplementations.SparseDotDenseVector(nnz, valuesX, indicesX, offsetX, vectorY, offsetY);
        #endregion

        #region Sparse BLAS Level 2

        public void Dcscgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] colOffsetsA, int[] rowIndicesA,
            double[] vectorX, int offsetX, double[] vectorY, int offsetY)
            => SparseBlasImplementations.CsrTimesVector(!transposeA, numColsA, numRowsA, valuesA, colOffsetsA, rowIndicesA,
                vectorX, offsetX, vectorY, offsetY);

        public void Dcsrgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA, int[] colIndicesA,
            double[] vectorX, int offsetX, double[] vectorY, int offsetY)
            => SparseBlasImplementations.CsrTimesVector(transposeA, numRowsA, numColsA, valuesA, rowOffsetsA, colIndicesA,
                vectorX, offsetX, vectorY, offsetY);

        /// <summary>
        /// Matrix vector multiplication y = A * x, with A being a symmetric matrix in Skyline format, where only the upper 
        /// triangle is stored.
        /// </summary>
        public void Dskymv(int order, double[] valuesA, int[] diagOffsetsA, double[] vectorX, double[] vectorY)
            => SparseBlasImplementations.SkylineTimesVector(order, valuesA, diagOffsetsA, vectorX, vectorY);

        /// <summary>
        /// Linear system solution x = inv(A) * b, with A being with a symmetric matrix in Skyline format, where only the upper 
        /// triangle is stored.
        /// </summary>
        public void Dskysv(int order, double[] valuesA, int[] diagOffsetsA, double[] vectorB, double[] vectorX)
            => SparseBlasImplementations.SkylineSystemSolution(order, valuesA, diagOffsetsA, vectorB, vectorX);

        /// <summary>
        /// Linear system solution x = inv(A) * b, with A being with a symmetric matrix in Skyline format, where only the upper 
        /// triangle is stored.
        /// </summary>
        public void Dskysv(int order, double[] valuesA, int[] diagOffsetsA, double[] vectorB, int offsetB,
            double[] vectorX, int offsetX)
            => SparseBlasImplementations.SkylineSystemSolutionWithOffsets(order, valuesA, diagOffsetsA,
                vectorB, offsetB, vectorX, offsetX);
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

            //TODO: This implementation uses level 2 BLAS. I should implement it from scratch as a lvl 3 BLAS, but is it worth 
            //      it without parallelism?
            for (int j = 0; j < numColsB; ++j)
            {
                Dcsrgemv(transposeA, numRowsA, numColsA, valuesA, rowOffsetsA, colIndicesA, matrixB, j * ldB, matrixC, j * ldC);
            }
        }

        public void Dskysm(int order, int numRhs, double[] valuesA, int[] diagOffsetsA, double[] vectorB, double[] vectorX)
        {
            //TODO: This implementation uses level 2 BLAS. I should implement it from scratch as a lvl 3 BLAS, but is it worth 
            //      it without parallelism?
            for (int j = 0; j < numRhs; ++j)
            {
                int offset = j * order;
                Dskysv(order, valuesA, diagOffsetsA, vectorB, offset, vectorX, offset);
            }
        }
        #endregion
    }
}
