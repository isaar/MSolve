using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    internal enum MatrixLayout
    {
        /// <summary>
        /// Row major layout.
        /// </summary>
        RowMajor,

        /// Column major layout.
        ColMajor
    }

    internal enum TransposeMatrix
    {
        /// <summary>
        /// The stored matrix will be accessed normally.
        /// </summary>
        NoTranspose,

        /// <summary>
        /// The transpose of the stored matrix will be accessed, without explicitly transposing.
        /// </summary>
        Transpose

        //Not supported yet.
        ///// <summary>
        ///// For complex matrices, the conjugate transpose of the stored matrix will be accessed, without explicitly transposing.
        ///// </summary>
        //ConjugateTranspose
    }

    internal enum StoredTriangle
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

    internal enum DiagonalValues
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

    internal enum MultiplicationSide
    {
        /// <summary>
        /// The corresponding matrix A will be multiplied on the left: A * whatever.
        /// </summary>
        Left,

        /// <summary>
        /// The corresponding matrix A will be multiplied on the right: whatever * A.
        /// </summary>
        Right
    }

    internal static class Translations
    {
        internal static string Translate(this MatrixLayout layout) => (layout == MatrixLayout.ColMajor) ? "F" : "C";
        internal static string Translate(this TransposeMatrix trans) => (trans == TransposeMatrix.Transpose) ? "T" : "N";
        internal static string Translate(this StoredTriangle uplo) => (uplo == StoredTriangle.Upper) ? "U" : "L";
        internal static string Translate(this DiagonalValues diag) => (diag == DiagonalValues.Unit) ? "U" : "N";
        internal static string Translate(this MultiplicationSide side) => (side == MultiplicationSide.Left) ? "L" : "R";
    }
}
