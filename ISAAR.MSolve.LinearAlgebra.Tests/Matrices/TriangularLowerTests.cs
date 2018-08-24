using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="TriangularLower"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class TriangularLowerTests
    {
        public const int order = 10;

        /// <summary>
        /// An invertible 10-by-10 lower triangular matrix.
        /// </summary>
        public static readonly double[,] invertible = new double[,] {
            {4.2289,         0,         0,         0,         0,         0,         0,         0,         0,         0},
            {0.9423,    6.5445,         0,         0,         0,         0,         0,         0,         0,         0},
            {5.9852,    4.0762,    0.9082,         0,         0,         0,         0,         0,         0,         0},
            {4.7092,    8.1998,    2.6647,    9.5769,         0,         0,         0,         0,         0,         0},
            {6.9595,    7.1836,    1.5366,    2.4071,    3.4446,         0,         0,         0,         0,         0},
            {6.9989,    9.6865,    2.8101,    6.7612,    7.8052,    7.7016,         0,         0,         0,         0},
            {6.3853,    5.3133,    4.4009,    2.8906,    6.7533,    3.2247,    1.9175,         0,         0,         0},
            {0.3360,    3.2515,    5.2714,    6.7181,    0.0672,    7.8474,    7.3843,    5.4659,         0,         0},
            {0.6881,    1.0563,    4.5742,    6.9514,    6.0217,    4.7136,    2.4285,    4.2573,    6.0730,         0},
            {3.1960,    6.1096,    8.7537,    0.6799,    3.8677,    0.3576,    9.1742,    6.4444,    4.5014,    6.1346}};

        /// <summary>
        /// A singular 10-by-10 lower triangular matrix.
        /// </summary>
        public static readonly double[,] singular = new double[,] {
            {4.2289,         0,         0,         0,         0,         0,         0,         0,         0,         0},
            {0.9423,    6.5445,         0,         0,         0,         0,         0,         0,         0,         0},
            {5.9852,    4.0762,    0.9082,         0,         0,         0,         0,         0,         0,         0},
            {4.7092,    8.1998,    2.6647,    9.5769,         0,         0,         0,         0,         0,         0},
            {6.9595,    7.1836,    1.5366,    2.4071,    3.4446,         0,         0,         0,         0,         0},
            {6.9595,    7.1836,    1.5366,    2.4071,    3.4446,         0,         0,         0,         0,         0},
            {6.3853,    5.3133,    4.4009,    2.8906,    6.7533,    3.2247,    1.9175,         0,         0,         0},
            {0.3360,    3.2515,    5.2714,    6.7181,    0.0672,    7.8474,    7.3843,    5.4659,         0,         0},
            {0.6881,    1.0563,    4.5742,    6.9514,    6.0217,    4.7136,    2.4285,    4.2573,    6.0730,         0},
            {3.1960,    6.1096,    8.7537,    0.6799,    3.8677,    0.3576,    9.1742,    6.4444,    4.5014,    6.1346}};

        /// <summary>
        /// A vector with length = 10.
        /// </summary>
        public static readonly double[] vector10 = 
            { 0.5822, 0.5407, 0.8699, 0.2648, 0.3181, 0.1192, 0.9398, 0.6456, 0.4795, 0.6393 };

        private static readonly Comparer comparer = new Comparer(1E-13);

        [Fact]
        public static void TestArrayCopy()
        {
            // invertible
            var A1 = TriangularLower.CreateFromArray(invertible);
            comparer.AssertEqual(invertible, A1.CopyToArray2D());

            // singular
            var A2 = TriangularLower.CreateFromArray(singular);
            comparer.AssertEqual(singular, A2.CopyToArray2D());
        }

        [Fact]
        public static void TestMatrixVectorMultiplication()
        {
            // invertible
            var A1 = TriangularLower.CreateFromArray(invertible);
            var x1 = Vector.CreateFromArray(vector10);
            var b1Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(invertible, vector10));
            Vector b1Computed = A1.MultiplyRight(x1);
            comparer.AssertEqual(b1Expected, b1Computed);

            // singular
            var A2 = TriangularLower.CreateFromArray(singular);
            var x2 = Vector.CreateFromArray(vector10);
            var b2Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(singular, vector10));
            Vector b2Computed = A2.MultiplyRight(x1);
            comparer.AssertEqual(b2Expected, b2Computed);
        }

        [Fact]
        public static void TestSystemSolution()
        {
            // invertible
            var A1 = TriangularLower.CreateFromArray(invertible);
            var b1 = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(invertible, vector10));
            var x1Expected = Vector.CreateFromArray(vector10);
            Vector x1Computed = A1.SolveLinearSystem(b1);
            comparer.AssertEqual(x1Expected, x1Computed);

            // singular
            var A2 = TriangularLower.CreateFromArray(singular);
            var b2 = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(singular, vector10));
            var x2Expected = Vector.CreateFromArray(vector10);
            Vector x2Computed = A2.SolveLinearSystem(b2);
            Assert.False(comparer.AreEqual(x2Expected, x2Computed));
        }

        [Fact]
        public static void TestTransposition()
        {
            // invertible
            var A1 = TriangularLower.CreateFromArray(invertible);
            var A1TransposeExpected = MatrixOperations.Transpose(invertible);
            TriangularUpper A1TransposeComputed = A1.Transpose(false);
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // singular
            var A2 = TriangularLower.CreateFromArray(singular);
            var A2TransposeExpected = MatrixOperations.Transpose(singular);
            TriangularUpper A2TransposeComputed = A2.Transpose(false);
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
