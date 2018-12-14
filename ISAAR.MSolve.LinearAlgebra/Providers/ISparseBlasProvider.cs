using System;
using System.Collections.Generic;
using System.Text;

//TODO: this and its implementations should be internal. The user should select a Provider that will then specify different BLAS,
//      SparseBLAS, LAPACK, etc. providers. The implementations should also be singletons or enums.
//TODO: some dimensions are redundant, since they can be read from the indexing arrays.
//TODO: Add offsets.
//TODO: For now this is my own interface. However there are some semi-standard ones. E.g. the inspector-executor interface, 
//      which also provides functions for A^T * B * A, a NIST interface, the default and simplified deprecated interfaces of MKL.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public interface ISparseBlasProvider
    {
        //TODO: also dot product
        void Daxpyi(int nnz, double[] alpha, double[] x, int offsetX, int[] indicesX, double[] y, int offsetY);

        /// <summary>
        /// Matrix-matrix multiplication C = op(A) * B, with A being a general matrix in 3-array CSC format (zero indexing) and
        /// B, C general matrices in full column major format.
        /// </summary>
        void Dcscgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] colOffsetsA, 
            int[] rowIndicesA, double[] b, double[] c);

        /// <summary>
        /// Matrix-vector multiplication y = A*x, with A being a general matrix in 3-array CSC format (zero indexing).
        /// </summary>
        void Dcscgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] colOffsetsA, int[] rowIndicesA, 
            double[] x, int offsetX, double[] y, int offsetY);

        /// <summary>
        /// Matrix-matrix multiplication C = op(A) * B, with A being a general matrix in 3-array CSR format (zero indexing) and
        /// B, C general matrices in full column major format.
        /// </summary>
        void Dcsrgemm(bool transposeA, int numRowsA, int numColsB, int numColsA, double[] valuesA, int[] rowOffsetsA,
            int[] colIndicesA, double[] b, double[] c);

        /// <summary>
        /// Matrix-vector multiplication y = A*x, with A being a general matrix in 3-array CSR format (zero indexing).
        /// </summary>
        void Dcsrgemv(bool transposeA, int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA, int[] colIndicesA,
            double[] x, int offsetX, double[] y, int offsetY);
    }
}
