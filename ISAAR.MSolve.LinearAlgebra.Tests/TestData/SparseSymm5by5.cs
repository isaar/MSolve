namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// A 5-by-5 sparse symmetric matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SparseSymm5by5
    {
        internal const int Order = 5;

        internal static double[,] Matrix => new double[,] {
            { 11.0,  0.0,  0.0,  0.0, 12.0 },
            {  0.0, 44.0, 45.0, 34.0,  0.0 },
            {  0.0, 45.0, 55.0,  0.0,  0.0 },
            {  0.0, 34.0,  0.0, 33.0, 23.0 },
            { 12.0,  0.0,  0.0, 23.0, 24.0 }
        };

        internal static double[] CscValues => new double[] {
            11.0, 12.0, 44.0, 45.0, 34.0, 45.0, 55.0, 34.0, 33.0, 23.0, 12.0, 23.0, 22.0 };
        internal static int[] CscRowIndices => new int[] { 0, 4, 1, 2, 3, 1, 2, 1, 3, 4, 0, 3, 4 };
        internal static int[] CscColOffsets => new int[] { 0, 2, 5, 7, 10, 13 };

        internal static double[] CscUpperValues => new double[] {
            11.0, 44.0, 45.0, 55.0, 34.0, 33.0, 12.0, 23.0, 22.0 };
        internal static int[] CscUpperRowIndices => new int[] { 0, 1, 1, 2, 1, 3, 0, 3, 4 };
        internal static int[] CscUpperColOffsets => new int[] { 0, 1, 2, 4, 6, 9 };

        /// <summary>
        /// New-to-old permutation resulting from Matlab's AMD algorithm. 0-based indexing is used.
        /// </summary>
        internal static int[] MatlabPermutationAMD => new int[] { 2, 1, 3, 0, 4 }; 
    }
}
