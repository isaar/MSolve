namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// A 5-by-5 sparse symmetric matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SparseSymm5by5
    {
        internal const int order = 5;

        internal static readonly double[,] matrix = {
            { 11.0,  0.0,  0.0,  0.0, 12.0 },
            {  0.0, 44.0, 45.0, 34.0,  0.0 },
            {  0.0, 45.0, 55.0,  0.0,  0.0 },
            {  0.0, 34.0,  0.0, 33.0, 23.0 },
            { 12.0,  0.0,  0.0, 23.0, 24.0 }
        };

        internal static readonly double[] cscValues = {
            11.0, 12.0, 44.0, 45.0, 34.0, 45.0, 55.0, 34.0, 33.0, 23.0, 12.0, 23.0, 22.0 };
        internal static readonly int[] cscRowIndices = { 0, 4, 1, 2, 3, 1, 2, 1, 3, 4, 0, 3, 4 };
        internal static readonly int[] cscColOffsets = { 0, 2, 5, 7, 10, 13 };

        internal static readonly double[] cscUpperValues = {
            11.0, 44.0, 45.0, 55.0, 34.0, 33.0, 12.0, 23.0, 22.0 };
        internal static readonly int[] cscUpperRowIndices = { 0, 1, 1, 2, 1, 3, 0, 3, 4 };
        internal static readonly int[] cscUpperColOffsets = { 0, 1, 2, 4, 6, 9 };

        /// <summary>
        /// New-to-old permutation resulting from Matlab's AMD algorithm. 0-based indexing is used.
        /// </summary>
        internal static readonly int[] matlabPermutationAMD = { 2, 1, 3, 0, 4 }; 
    }
}
