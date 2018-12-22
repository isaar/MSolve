using System.Runtime.CompilerServices;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Describes how a matrix (2D) is stored in the system's memory as a 1D array.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal enum MatrixLayout
    {
        /// <summary>
        /// Row major layout. Entries corresponding to consecutive columns are stored consecutively. This is the default in C 
        /// and C#. 
        /// </summary>
        RowMajor,

        /// Column major layout. Entries corresponding to consecutive rows are stored consecutively. This is the default in 
        /// FORTRAN and most linear algebra libraries that use the BLAS & LAPACK interfaces.
        ColMajor
    }

    /// <summary>
    /// Describes whether a matrix will be transposed or not, during a BLAS or LAPACK operation. Note that if 
    /// <see cref="TransposeMatrix.Transpose"/> is selected, they matrix is not explicitly transposed (usually). Instead access
    /// to the matrix entries is transposed, e.g. A[j, i] instead of A[i, j].
    /// Authors: Serafeim Bakalakos
    /// </summary>
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

        //Not supported yet. When it is supported, I must check all code where the exact value of TransposeMatrix is searched.
        ///// <summary>
        ///// For complex matrices, the conjugate transpose of the stored matrix will be accessed, without explicitly transposing.
        ///// </summary>
        //ConjugateTranspose
    }

    /// <summary>
    /// Describes which triangle of a triangular or symmetric matrix is explicitly stored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal enum StoredTriangle
    {
        /// <summary>
        /// The upper triangle of the matrix is stored, including the diagonal.
        /// </summary>
        Upper,

        /// <summary>
        /// The lower triangle of the matrix is stored, including the diagonal.
        /// </summary>
        Lower
    }

    /// <summary>
    /// Describes whether the diagonal values of a matrix are explicitly stored or not. Applicable to triangular matrices only.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal enum DiagonalValues
    {
        /// <summary>
        /// The entries of the diagonal are assumed to be 1. 
        /// </summary>
        Unit,

        /// <summary>
        /// The entries of the diagonal are arbitrary numbers explicitly stored. 
        /// </summary>
        NonUnit
    }

    /// <summary>
    /// Describes on which side of a matrix-matrix multiplication is the matrix of interest located.
    /// Authors: Serafeim Bakalakos
    /// </summary>
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

    /// <summary>
    /// Extension methods to convert BLAS & LAPACK oriented enums into the strings expected by the actual BLAS & LAPACK 
    /// methods.
    /// Authors: Serafeim Bakalakos
    /// </summary>
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
