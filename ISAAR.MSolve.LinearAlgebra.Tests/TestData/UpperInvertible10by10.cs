using System.IO;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// An invertible 10-by-10 upper triangular matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class UpperInvertible10by10
    {
        internal const int order = 10;

        internal static readonly double[,] matrix = new double[,] {
            {4.2289,    5.3086,   7.7880,    5.1805,    2.5479,    9.1599,    1.7587,    2.6906,    6.4762,    4.5873},
            { 0,        6.5445,   4.2345,    9.4362,    2.2404,    0.0115,    7.2176,    7.6550,    6.7902,    6.6194},
            { 0,         0,       0.9082,    6.3771,    6.6783,    4.6245,    4.7349,    1.8866,    6.3579,    7.7029},
            { 0,         0,        0,        9.5769,    8.4439,    4.2435,    1.5272,    2.8750,    9.4517,    3.5022},
            { 0,         0,        0,        0,         3.4446,    4.6092,    3.4112,    0.9111,    2.0893,    6.6201},
            { 0,         0,        0,        0,         0,         7.7016,    6.0739,    5.7621,    7.0928,    4.1616},
            { 0,         0,        0,        0,         0,          0,        1.9175,    6.8336,    2.3623,    8.4193},
            { 0,         0,        0,        0,         0,          0,        0,         5.4659,    1.1940,    8.3292},
            { 0,         0,        0,        0,         0,          0,        0,         0,         6.0730,    2.5644},
            { 0,         0,        0,        0,         0,          0,        0,         0,         0,         6.1346}
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
            + @"\Resources\UpperInvertible10by10.txt";
    }
}
