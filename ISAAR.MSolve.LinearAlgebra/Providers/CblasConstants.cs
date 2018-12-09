using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public enum CblasLayout
    {
        /// <summary>
        /// Row major layout.
        /// </summary>
        RowMajor,

        /// Column major layout.
        ColMajor
    }

    public enum CblasTranspose
    {
        /// <summary>
        /// The stored matrix will be accessed normally.
        /// </summary>
        NoTranspose,

        /// <summary>
        /// The transpose of the stored matrix will be accessed, without explicitly transposing.
        /// </summary>
        Transpose,

        /// <summary>
        /// For complex matrices, the conjugate transpose of the stored matrix will be accessed, without explicitly transposing.
        /// </summary>
        ConjugateTranspose
    }

    public enum CblasTriangular
    {
        /// <summary>
        /// The stored matrix is upper triangular.
        /// </summary>
        Upper,

        /// <summary>
        /// The stored matrix is lower triangular.
        /// </summary>
        Lower
    }

    public enum CblasDiagonal
    {
        /// <summary>
        /// The entries of the diagonal are 1. Applicable to triangular matrices only. Therefore the matrix is upper unit 
        /// triangular or lower unit triangular.
        /// </summary>
        Unit,

        /// <summary>
        /// The entries of the diagonal are arbitrary numbers. Applicable to triangular matrices only.
        /// </summary>
        NonUnit
    }
}
