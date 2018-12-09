using System.Runtime.CompilerServices;
using IntelMKL.LP64;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    internal class MklConstants
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static CBLAS_LAYOUT Translate(CblasLayout layout)
            => layout == CblasLayout.ColMajor ? CBLAS_LAYOUT.CblasColMajor : CBLAS_LAYOUT.CblasRowMajor;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static CBLAS_TRANSPOSE Translate(CblasTranspose trans)
        {
            if (trans == CblasTranspose.NoTranspose) return CBLAS_TRANSPOSE.CblasNoTrans;
            else if (trans == CblasTranspose.Transpose) return CBLAS_TRANSPOSE.CblasTrans;
            else return CBLAS_TRANSPOSE.CblasConjTrans;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static CBLAS_UPLO Translate(CblasTriangular uplo)
            => uplo == CblasTriangular.Upper ? CBLAS_UPLO.CblasUpper : CBLAS_UPLO.CblasLower;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static CBLAS_DIAG Translate(CblasDiagonal diag)
            => diag == CblasDiagonal.Unit ? CBLAS_DIAG.CblasUnit : CBLAS_DIAG.CblasNonUnit;
    }
}
