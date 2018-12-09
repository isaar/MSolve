using System;
using System.Runtime.CompilerServices;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    internal class ManagedConstants
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(CBlasLayout layout, CBlasTranspose trans) //TODO: Not sure if this covers every case.
        {
            if (layout == CBlasLayout.ColMajor)
            {
                if (trans == CBlasTranspose.NoTranspose) return "N";
                else if (trans == CBlasTranspose.Transpose) return "T";
                else return "C";
            }
            else // BLAS does not have a row major layout, so we need to transpose
            {
                if (trans == CBlasTranspose.NoTranspose) return "T";
                else if (trans == CBlasTranspose.Transpose) return "N";
                else throw new ArgumentException(
                    "The combination of row major layout and conjugate transpose access is not supported");
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(CBlasTriangular uplo)
            => uplo == CBlasTriangular.Upper ? "U" : "L";

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(CBlasDiagonal diag)
            => diag == CBlasDiagonal.Unit ? "U" : "N";
    }
}
