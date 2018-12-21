using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;

//TODO: these should be internal
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public enum MatrixLayout
    {
        /// <summary>
        /// Row major layout.
        /// </summary>
        RowMajor,

        /// Column major layout.
        ColMajor
    }

    public enum TransposeMatrix
    {
        /// <summary>
        /// The stored matrix will be accessed normally.
        /// </summary>
        NoTranspose,

        /// <summary>
        /// The transpose of the stored matrix will be accessed, without explicitly transposing.
        /// </summary>
        Transpose

        //Not supported yet. When it is supported, I must check all code where the exact value of TransposeMatrix is searched.
        ///// <summary>
        ///// For complex matrices, the conjugate transpose of the stored matrix will be accessed, without explicitly transposing.
        ///// </summary>
        //ConjugateTranspose
    }

    public enum StoredTriangle
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

    public enum DiagonalValues
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

    public enum MultiplicationSide
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
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(this MatrixLayout layout) => (layout == MatrixLayout.ColMajor) ? "F" : "C";

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(this TransposeMatrix trans) => (trans == TransposeMatrix.Transpose) ? "T" : "N";

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(this StoredTriangle uplo) => (uplo == StoredTriangle.Upper) ? "U" : "L";

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(this DiagonalValues diag) => (diag == DiagonalValues.Unit) ? "U" : "N";

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(this MultiplicationSide side) => (side == MultiplicationSide.Left) ? "L" : "R";
    }
}
