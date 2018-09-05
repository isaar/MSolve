using System.IO;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// A 10-by-5 sparse matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SparseRectangular10by5
    {
        internal static readonly int numRows = 10;
        internal static readonly int numCols = 5;

        internal static readonly double[,] matrix = {
            { 0.00000,   0.00000,   1.20805,   2.09065,   0.00000 },
            { 0.72602,   3.72194,   0.64231,   2.74205,   0.00000 },
            { 1.68190,   0.00000,   2.78180,   0.00000,   0.02756 },
            { 0.00000,   0.00000,   0.83011,   0.00000,   0.54185 },
            { 2.44372,   1.52758,   0.00000,   2.98743,   0.00000 },
            { 0.00000,   0.00000,   0.77108,   3.53890,   0.00000 },
            { 0.53085,   0.00000,   0.00000,   3.04244,   0.00000 },
            { 0.00000,   0.00000,   0.00000,   0.00000,   3.85200 },
            { 0.42069,   0.00000,   2.65133,   2.70976,   0.52623 },
            { 3.19534,   0.00000,   0.00000,   3.80848,   0.00000 }
        };

        internal static readonly double[] csrValues = {
            1.20805, 2.09065,
            0.72602, 3.72194, 0.64231, 2.74205,
            1.68190, 2.78180, 0.02756,
            0.83011, 0.54185,
            2.44372, 1.52758, 2.98743,
            0.77108, 3.53890,
            0.53085, 3.04244,
            3.85200,
            0.42069, 2.65133, 2.70976, 0.52623,
            3.19534, 3.80848
        };

        internal static readonly int[] csrColIndices = {
            2, 3,
            0, 1, 2, 3,
            0, 2, 4,
            2, 4,
            0, 1, 3,
            2, 3,
            0, 3,
            4,
            0, 2, 3, 4,
            0, 3
        };

        internal static readonly int[] csrRowOffsets = { 0, 2, 6, 9, 11, 14, 16, 18, 19, 23, 25 };


        internal static readonly double[] cscValues = {
            0.72602, 1.68190, 2.44372, 0.53085, 0.42069, 3.19534,
            3.72194, 1.52758,
            1.20805, 0.64231, 2.78180, 0.83011, 0.77108, 2.65133,
            2.09065, 2.74205, 2.98743, 3.53890, 3.04244, 2.70976, 3.80848,
            0.02756, 0.54185, 3.85200, 0.52623
        };

        internal static readonly int[] cscRowIndices = {
            1, 2, 4, 6, 8, 9,
            1, 4,
            0, 1, 2, 3, 5, 8,
            0, 1, 4, 5, 6, 8, 9,
            2, 3, 7, 8
        };

        internal static readonly int[] cscColOffsets = { 0, 6, 8, 14, 21, 25 };

        /// <summary>
        /// A vector with length = 5, such that <see cref="matrix"/> * <see cref="lhs5"/> = <see cref="rhs10"/>.
        /// </summary>
        internal static readonly double[] lhs5 = { 7.5025, 9.8100, 2.3352, 0.9623, 3.8458 };

        /// <summary>
        /// A vector with length = 10, such that <see cref="matrix"/> * <see cref="lhs5"/> = <see cref="rhs10"/>.
        /// </summary>
        internal static readonly double[] rhs10 = {
            4.832870855, 46.097793477, 19.220504358, 4.022319602, 36.194372989, 5.206109486, 6.910442137, 14.814021600, 13.978989923, 27.637938654 };

        /// <summary>
        /// A vector with length = 10, such that <see cref="matrix"/> * <see cref="lhs10"/> = <see cref="rhs5"/>.
        /// </summary>
        internal static readonly double[] lhs10 = {
            7.5025, 9.8100, 2.3352, 0.9623, 3.8458, 5.0027, 5.7026, 9.7663, 4.9286, 4.0088 };

        /// <summary>
        /// A vector with length = 5, such that <see cref="matrix"/> * <see cref="lhs10"/> = <see cref="rhs5"/>.
        /// </summary>
        internal static readonly double[] rhs5 = { 38.358004392, 42.386998564, 39.584157392, 117.750301553, 40.799145145 };

        internal static readonly string fullFormatPath =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparseRectangular10by5_FullFormat.txt";

        internal static readonly string csrArraysPath =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparseRectangular10by5_CsrArrays.txt";

        internal static readonly string cscArraysPath =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparseRectangular10by5_CscArrays.txt";

        /// <summary>
        /// The entries were iteration in CSC order.
        /// </summary>
        internal static readonly string coordinateFormatPath =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparseRectangular10by5_CoordinateFormat.txt";
    }
}
