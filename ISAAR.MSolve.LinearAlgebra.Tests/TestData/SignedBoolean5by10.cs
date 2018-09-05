using System.IO;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// A 5-by-10 sparse matrix, that only contains 0, 1, -1.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SignedBoolean5by10
    {
        internal static readonly int numRows = 5;
        internal static readonly int numCols = 10;

        internal static readonly double[,] A1 =
        {
            { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0,-1, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0,-1, 0, 0 }
        };

        internal static readonly double[,] A2 =
        {
            { 1,-1, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 1,-1, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 1,-1, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 1,-1, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 1,-1, 0, 0, 0, 0 }
        };

        internal static readonly double[] x10 = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.1 };
        internal static readonly double[] A1TimesX10 = { 5.5, 0.0, -2.2, 0.0, -8.8 };
        internal static readonly double[] A2TimesX10 = { -1.1, -1.1, -1.1, -1.1, -1.1 };
        internal static readonly double[] x5 = { 1.1, 2.2, 3.3, 4.4, 5.5 };
        internal static readonly double[] transpA1TimesX5 = { 0.0, -3.3, 0.0, 0.0, 1.1, 0.0, 0.0, -5.5, 0.0, 0.0 };
        internal static readonly double[] transpA2TimesX5 = { 1.1, 1.1, 1.1, 1.1, 1.1, -5.5, 0.0, 0.0, 0.0, 0.0, };

        internal static readonly string filePath1 = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
           + @"\Resources\SignedBoolean1.txt";

        internal static readonly string filePath2 = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
           + @"\Resources\SignedBoolean2.txt";
    }
}
