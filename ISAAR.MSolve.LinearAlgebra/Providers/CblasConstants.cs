using System;
using System.Collections.Generic;
using System.Text;

//TODO: these and the Translate methods shoudl be refactored. Also note that I will only use BLAS, LAPACK (for starters), 
//      which means that the translate methods are the same and thus can be written here using enum classes.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public enum CBlasLayout
    {
        /// <summary>
        /// Row major layout.
        /// </summary>
        RowMajor,

        /// Column major layout.
        ColMajor
    }

    public enum CBlasTranspose
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

    public enum CBlasTriangular
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

    public enum CBlasDiagonal
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
