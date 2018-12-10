using System.IO;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// A symmetric 10-by-10 matrix that is singular.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SymmSingular10by10
    {
        internal const int Order = 10;

        internal static double[,] Matrix => new double[,] {
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300}};

        internal static double[] Lhs => new double[] {
            2.6621, 3.5825, 0.8965, 1.6827, 0.9386, 1.6096, 2.0193, 2.7428, 0.2437, 2.7637
        };

        internal static double[] Rhs => MatrixOperations.MatrixTimesVector(Matrix, Lhs);

        internal static string FilePath => Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SymmSingular10by10.txt";
    }
}
