using System;
using System.Runtime.CompilerServices;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    internal class ManagedConstants
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(CblasLayout layout, CblasTranspose trans) //TODO: Not sure if this covers every case.
        {
            if (layout == CblasLayout.ColMajor)
            {
                if (trans == CblasTranspose.NoTranspose) return "N";
                else if (trans == CblasTranspose.Transpose) return "T";
                else return "C";
            }
            else // BLAS does not have a row major layout, so we need to transpose
            {
                if (trans == CblasTranspose.NoTranspose) return "T";
                else if (trans == CblasTranspose.Transpose) return "N";
                else throw new ArgumentException(
                    "The combination of row major layout and conjugate transpose access is not supported");
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(CblasTriangular uplo)
            => uplo == CblasTriangular.Upper ? "U" : "L";

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static string Translate(CblasDiagonal diag)
            => diag == CblasDiagonal.Unit ? "U" : "N";
    }
}
