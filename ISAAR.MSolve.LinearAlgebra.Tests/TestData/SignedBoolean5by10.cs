using System.IO;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// A 5-by-10 sparse matrix, that only contains 0, 1, -1.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SignedBoolean5by10
    {
        internal const int NumRows = 5;
        internal const int NumCols = 10;

        internal static double[,] A1 => new double[,]
        {
            { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0,-1, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0,-1, 0, 0 }
        };

        internal static double[,] A2 => new double[,]
        {
            { 1,-1, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 1,-1, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 1,-1, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 1,-1, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 1,-1, 0, 0, 0, 0 }
        };

        internal static double[] X10 => new double[] { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.1 };
        internal static double[] A1TimesX10 => new double[] { 5.5, 0.0, -2.2, 0.0, -8.8 };
        internal static double[] A2TimesX10 => new double[] { -1.1, -1.1, -1.1, -1.1, -1.1 };
        internal static double[] X5 => new double[] { 1.1, 2.2, 3.3, 4.4, 5.5 };
        internal static double[] TranspA1TimesX5 => new double[] { 0.0, -3.3, 0.0, 0.0, 1.1, 0.0, 0.0, -5.5, 0.0, 0.0 };
        internal static double[] TranspA2TimesX5 => new double[] { 1.1, 1.1, 1.1, 1.1, 1.1, -5.5, 0.0, 0.0, 0.0, 0.0 };

        internal static string FilePath1 => Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
           + @"\Resources\SignedBoolean1.txt";

        internal static string FilePath2 => Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
           + @"\Resources\SignedBoolean2.txt";
    }
}
