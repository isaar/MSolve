//TODO: some dimensions are redundant, since they can be read from the indexing arrays.
//TODO: Add offsets.
//TODO: For now this is my own interface. However there are some semi-standard ones. E.g. the inspector-executor interface, 
//      which also provides functions for A^T * B * A, a NIST interface, the default and simplified deprecated interfaces of MKL.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Provides linear algebra operations for sparse matrices and vectors, similar to the ones defined by BLAS for dense 
    /// matrices and vectors (see <see cref="IBlasProvider"/>).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal interface ISparseBlasProvider
    {
        #region Sparse BLAS Level 1

        /// <summary>
        /// Linear combination y = alpha * x + y, where x is a sparse vector and y is a dense vector.
        /// </summary>
        void Daxpyi(int nnz, double alpha, double[] valuesX, int[] indicesX, int offsetX, double[] vectorY, int offsetY);

        /// <summary>
        /// Dot product x * y, where x is a sparse vector and y is a dense vector.
        /// </summary>
        double Ddoti(int nnz, double[] valuesX, int[] indicesX, int offsetX, double[] vectorY, int offsetY);
        #endregion

        #region Sparse BLAS Level 2

        /// <summary>
        /// Matrix-vector multiplication y = A*x, with A being a general matrix in 3-array CSC format (zero indexing).
        /// </summary>
        void Dcscgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] colOffsetsA, int[] rowIndicesA, 
            double[] vectorX, int offsetX, double[] vectorY, int offsetY);

        /// <summary>
        /// Matrix-vector multiplication y = A*x, with A being a general matrix in 3-array CSR format (zero indexing).
        /// </summary>
        void Dcsrgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA, int[] colIndicesA,
            double[] vectorX, int offsetX, double[] vectorY, int offsetY);
        #endregion

        #region Sparse BLAS Level 3

        /// <summary>
        /// Matrix-matrix multiplication C = op(A) * B, with A being a general matrix in 3-array CSC format (zero indexing) and
        /// B, C general matrices in full column major format.
        /// </summary>
        void Dcscgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] colOffsetsA,
            int[] rowIndicesA, double[] matrixB, double[] matrixC);

        /// <summary>
        /// Matrix-matrix multiplication C = op(A) * B, with A being a general matrix in 3-array CSR format (zero indexing) and
        /// B, C general matrices in full column major format.
        /// </summary>
        void Dcsrgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] rowOffsetsA,
            int[] colIndicesA, double[] matrixB, double[] matrixC);
        #endregion
    }
}
