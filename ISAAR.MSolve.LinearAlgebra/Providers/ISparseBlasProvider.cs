using System;
using System.Collections.Generic;
using System.Text;

//TODO: this and its implementations should be internal. The user should select a Provider that will then specify different BLAS,
//      SparseBLAS, LAPACK, etc. providers. The implementations shoud lalso be singletons or enums.
//TODO: some dimensions are redundant, since they can be read from the indexing arrays.
//TODO: Add offsets.
//TODO: For now this is my own interface. However there are some semi-standard ones. E.g. the inspector-executor interface, 
//      which also provides functions for A^T * B * A, a NIST interface, the default and simplified deprecated interfaces of MKL.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public interface ISparseBlasProvider
    {
        void Daxpyi(int nnz, double[] alpha, double[] x, int offsetX, int[] indicesX, double[] y, int offsetY);

        /// <summary>
        /// Matrix vector multiplication, with a general matrix in 3-array CSC format (zero indexing). See
        /// </summary>
        void Dcscgemv(bool transpose, int numRows, int numCols, double[] values, int[] colOffsets, int[] rowIndices, 
            double[] lhs, double[] rhs);

        /// <summary>
        /// Matrix vector multiplication, with a general matrix in 3-array CSR format (zero indexing). See
        /// </summary>
        void Dcsrgemv(bool transpose, int numRows, int numColumns, double[] values, int[] rowOffsets, int[] colIndices, 
            double[] lhs, double[] rhs);
    }
}
