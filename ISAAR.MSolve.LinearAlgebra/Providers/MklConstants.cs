using System.Runtime.CompilerServices;
using IntelMKL.LP64;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    internal class MklConstants
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static CBLAS_LAYOUT Translate(CBlasLayout layout)
            => layout == CBlasLayout.ColMajor ? CBLAS_LAYOUT.CblasColMajor : CBLAS_LAYOUT.CblasRowMajor;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static CBLAS_TRANSPOSE Translate(CBlasTranspose trans)
        {
            if (trans == CBlasTranspose.NoTranspose) return CBLAS_TRANSPOSE.CblasNoTrans;
            else if (trans == CBlasTranspose.Transpose) return CBLAS_TRANSPOSE.CblasTrans;
            else return CBLAS_TRANSPOSE.CblasConjTrans;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static CBLAS_UPLO Translate(CBlasTriangular uplo)
            => uplo == CBlasTriangular.Upper ? CBLAS_UPLO.CblasUpper : CBLAS_UPLO.CblasLower;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static CBLAS_DIAG Translate(CBlasDiagonal diag)
            => diag == CBlasDiagonal.Unit ? CBLAS_DIAG.CblasUnit : CBLAS_DIAG.CblasNonUnit;
    }
}
