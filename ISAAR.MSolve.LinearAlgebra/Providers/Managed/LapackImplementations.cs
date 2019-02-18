using System;

//TODO: Some of these could be done by calling other BLAS, LAPACK functions. See the LAPACK source.
//TODO: Use a port of reference LAPACK here. Custom implementations should be in another namespace (e.g. Factorizations or 
//      Triangulations)
namespace ISAAR.MSolve.LinearAlgebra.Providers.Managed
{
    /// <summary>
    /// Custom and unoptimized managed implementations of LAPACK like operations, for which I have not found 3rd party 
    /// implementations yet. If performance is of concern, then the user should choose the native optimized providers (e.g. MKL).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class LapackImplementations
    {
        internal static void CholeskyLowerFullColMajor(int n, double[] a, int offsetA, int ldA, ref int info)
            => CholeskyTemplate(n, a, ref info, (i, j) => offsetA + j * ldA + i);

        internal static void CholeskyLowerPackedColMajor(int n, double[] a, int offsetA, ref int info)
            => CholeskyTemplate(n, a, ref info, (i, j) => offsetA + i + (j * (2 * n + 1 - j)) / 2);

        internal static void CholeskyUpperFullColMajor(int n, double[] a, int offsetA, int ldA, ref int info)
            => CholeskyTemplate(n, a, ref info, (i, j) => offsetA + i * ldA + j); // Do not transpose access

        internal static void CholeskyUpperPackedColMajor(int n, double[] a, int offsetA, ref int info)
            => CholeskyTemplate(n, a, ref info, (i, j) => offsetA + j + (i * (i + 1)) / 2); //Transpose access

        /// <summary>
        /// For Cholesky algorithm, see 
        /// https://algowiki-project.org/en/Cholesky_decomposition#Software_implementation_of_the_algorithm. 
        /// For the indexing formulas see
        /// https://software.intel.com/en-us/mkl-developer-reference-c-matrix-storage-schemes-for-lapack-routines.
        /// This algorithm operates on L only. The caller might need to transpose (i, j), in order to operate on U instead.
        /// </summary>
        private static void CholeskyTemplate(int n, double[] a, ref int info, Func<int, int, int> FindIndex)
        {
            for (int i = 0; i < n; ++i)
            {
                // Calculate the diagonal entry of column i
                int diagIdx = FindIndex(i, i);
                double diagValue = a[diagIdx];
                for (int k = 0; k < i; ++k)
                {
                    int colIdx = FindIndex(i, k);
                    diagValue -= a[colIdx] * a[colIdx];
                }

                if (diagValue <= 0)
                {
                    info = i;
                    return;
                }
                diagValue = Math.Sqrt(diagValue);
                a[diagIdx] = diagValue;

                // Calculate the subdiagonal entries of column i
                for (int j = i+1; j < n; ++j)
                {
                    int rowIdx = FindIndex(j, i);
                    double nominator = a[rowIdx];
                    for (int k = 0; k < i; ++k)
                    {
                        nominator -= a[FindIndex(i, k)] * a[FindIndex(j, k)];
                    }
                    a[rowIdx] = nominator / diagValue;
                }
            }
            info = 0;
        }
    }
}
