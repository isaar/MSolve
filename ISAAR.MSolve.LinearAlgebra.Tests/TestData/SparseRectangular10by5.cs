using System.IO;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// A 10-by-5 sparse matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SparseRectangular10by5
    {
        internal const int NumRows = 10;
        internal const int NumCols = 5;

        internal static double[,] Matrix => new double[,] {
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

        internal static double[] CsrValues => new double[] {
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

        internal static int[] CsrColIndices => new int[] {
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

        internal static int[] CsrRowOffsets => new int[] { 0, 2, 6, 9, 11, 14, 16, 18, 19, 23, 25 };


        internal static double[] CscValues => new double[] {
            0.72602, 1.68190, 2.44372, 0.53085, 0.42069, 3.19534,
            3.72194, 1.52758,
            1.20805, 0.64231, 2.78180, 0.83011, 0.77108, 2.65133,
            2.09065, 2.74205, 2.98743, 3.53890, 3.04244, 2.70976, 3.80848,
            0.02756, 0.54185, 3.85200, 0.52623
        };

        internal static int[] CscRowIndices => new int[] {
            1, 2, 4, 6, 8, 9,
            1, 4,
            0, 1, 2, 3, 5, 8,
            0, 1, 4, 5, 6, 8, 9,
            2, 3, 7, 8
        };

        internal static int[] CscColOffsets => new int[] { 0, 6, 8, 14, 21, 25 };

        /// <summary>
        /// A vector with length = 5, such that <see cref="Matrix"/> * <see cref="Lhs5"/> = <see cref="Rhs10"/>.
        /// </summary>
        internal static double[] Lhs5 => new double[] { 7.5025, 9.8100, 2.3352, 0.9623, 3.8458 };

        /// <summary>
        /// A vector with length = 10, such that <see cref="Matrix"/> * <see cref="Lhs5"/> = <see cref="Rhs10"/>.
        /// </summary>
        internal static double[] Rhs10 => new double[] {
            4.832870855, 46.097793477, 19.220504358, 4.022319602, 36.194372989, 5.206109486, 6.910442137, 14.814021600, 13.978989923, 27.637938654 };

        /// <summary>
        /// A vector with length = 10, such that <see cref="Matrix"/> * <see cref="Lhs10"/> = <see cref="Rhs5"/>.
        /// </summary>
        internal static double[] Lhs10 => new double[] {
            7.5025, 9.8100, 2.3352, 0.9623, 3.8458, 5.0027, 5.7026, 9.7663, 4.9286, 4.0088 };

        /// <summary>
        /// A vector with length = 5, such that <see cref="Matrix"/> * <see cref="Lhs10"/> = <see cref="Rhs5"/>.
        /// </summary>
        internal static double[] Rhs5 => new double[] { 38.358004392, 42.386998564, 39.584157392, 117.750301553, 40.799145145 };

        internal static string FullFormatPath =>
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparseRectangular10by5_FullFormat.txt";

        internal static string CsrArraysPath =>
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparseRectangular10by5_CsrArrays.txt";

        internal static string CscArraysPath =>
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparseRectangular10by5_CscArrays.txt";

        /// <summary>
        /// The entries were iterated in CSC order.
        /// </summary>
        internal static string CoordinateFormatPath =>
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparseRectangular10by5_CoordinateFormat.txt";
    }
}
