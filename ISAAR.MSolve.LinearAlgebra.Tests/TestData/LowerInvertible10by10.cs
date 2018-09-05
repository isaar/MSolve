using System.IO;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// An invertible 10-by-10 lower triangular matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class LowerInvertible10by10
    {
        internal const int order = 10; 

        internal static readonly double[,] matrix = new double[,] {
            {4.2289,         0,         0,         0,         0,         0,         0,         0,         0,         0},
            {0.9423,    6.5445,         0,         0,         0,         0,         0,         0,         0,         0},
            {5.9852,    4.0762,    0.9082,         0,         0,         0,         0,         0,         0,         0},
            {4.7092,    8.1998,    2.6647,    9.5769,         0,         0,         0,         0,         0,         0},
            {6.9595,    7.1836,    1.5366,    2.4071,    3.4446,         0,         0,         0,         0,         0},
            {6.9989,    9.6865,    2.8101,    6.7612,    7.8052,    7.7016,         0,         0,         0,         0},
            {6.3853,    5.3133,    4.4009,    2.8906,    6.7533,    3.2247,    1.9175,         0,         0,         0},
            {0.3360,    3.2515,    5.2714,    6.7181,    0.0672,    7.8474,    7.3843,    5.4659,         0,         0},
            {0.6881,    1.0563,    4.5742,    6.9514,    6.0217,    4.7136,    2.4285,    4.2573,    6.0730,         0},
            {3.1960,    6.1096,    8.7537,    0.6799,    3.8677,    0.3576,    9.1742,    6.4444,    4.5014,    6.1346}
        };

        /// <summary>
        /// A vector with length = 10, such that <see cref="matrix"/> * <see cref="lhs"/> = <see cref="rhs"/>
        /// </summary>
        internal static readonly double[] lhs = {
            0.5822, 0.5407, 0.8699, 0.2648, 0.3181, 0.1192, 0.9398, 0.6456, 0.4795, 0.6393
        };

        /// <summary>
        /// A vector with length = 10, such that <see cref="matrix"/> * <see cref="lhs"/> = <see cref="rhs"/>
        /// </summary>
        internal static readonly double[] rhs = MatrixOperations.MatrixTimesVector(matrix, lhs);

        internal static readonly string filePath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName 
            + @"\Resources\LowerInvertible10by10.txt";
    }
}
