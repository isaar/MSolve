using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="TriangularUpper"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class TriangularUpperTests
    {
        public const int order = 10;

        /// <summary>
        /// An invertible 10-by-10 upper triangular matrix.
        /// </summary>
        public static readonly double[,] invertible = new double[,] {
            {4.2289,    5.3086,   7.7880,    5.1805,    2.5479,    9.1599,    1.7587,    2.6906,    6.4762,    4.5873},
            { 0,        6.5445,   4.2345,    9.4362,    2.2404,    0.0115,    7.2176,    7.6550,    6.7902,    6.6194},
            { 0,         0,       0.9082,    6.3771,    6.6783,    4.6245,    4.7349,    1.8866,    6.3579,    7.7029},
            { 0,         0,        0,        9.5769,    8.4439,    4.2435,    1.5272,    2.8750,    9.4517,    3.5022},
            { 0,         0,        0,        0,         3.4446,    4.6092,    3.4112,    0.9111,    2.0893,    6.6201},
            { 0,         0,        0,        0,         0,         7.7016,    6.0739,    5.7621,    7.0928,    4.1616},
            { 0,         0,        0,        0,         0,          0,        1.9175,    6.8336,    2.3623,    8.4193},
            { 0,         0,        0,        0,         0,          0,        0,         5.4659,    1.1940,    8.3292},
            { 0,         0,        0,        0,         0,          0,        0,         0,         6.0730,    2.5644},
            { 0,         0,        0,        0,         0,          0,        0,         0,         0,         6.1346}};
        
        /// <summary>
        /// A singular 10-by-10 upper triangular matrix.
        /// </summary>
        public static readonly double[,] singular = new double[,] {
            {4.2289,    5.3086,   7.7880,    5.1805,    2.5479,    9.1599,    1.7587,    2.6906,    6.4762,    4.5873},
            { 0,        6.5445,   4.2345,    9.4362,    2.2404,    0.0115,    7.2176,    7.6550,    6.7902,    6.6194},
            { 0,         0,       0.0000,    6.3771,    6.6783,    4.6245,    4.7349,    1.8866,    6.3579,    7.7029},
            { 0,         0,        0,        9.5769,    8.4439,    4.2435,    1.5272,    2.8750,    9.4517,    3.5022},
            { 0,         0,        0,        0,         3.4446,    4.6092,    3.4112,    0.9111,    2.0893,    6.6201},
            { 0,         0,        0,        0,         0,         7.7016,    6.0739,    5.7621,    7.0928,    4.1616},
            { 0,         0,        0,        0,         0,         0,         1.9175,    6.8336,    2.3623,    8.4193},
            { 0,         0,        0,        0,         0,         0,         0,         5.4659,    1.1940,    8.3292},
            { 0,         0,        0,        0,         0,         0,         0,         0,         6.0730,    2.5644},
            { 0,         0,        0,        0,         0,         0,         0,         0,         0,         6.1346}};

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
            var A1 = TriangularUpper.CreateFromArray(invertible);
            comparer.AssertEqual(invertible, A1.CopyToArray2D());

            // singular
            var A2 = TriangularUpper.CreateFromArray(singular);
            comparer.AssertEqual(singular, A2.CopyToArray2D());
        }

        [Fact]
        public static void TestMatrixVectorMultiplication()
        {
            // invertible
            var A1 = TriangularUpper.CreateFromArray(invertible);
            var x1 = Vector.CreateFromArray(vector10);
            var b1Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(invertible, vector10));
            Vector b1Computed = A1.MultiplyRight(x1);
            comparer.AssertEqual(b1Expected, b1Computed);

            // singular
            var A2 = TriangularUpper.CreateFromArray(singular);
            var x2 = Vector.CreateFromArray(vector10);
            var b2Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(singular, vector10));
            Vector b2Computed = A2.MultiplyRight(x1);
            comparer.AssertEqual(b2Expected, b2Computed);
        }

        [Fact]
        public static void TestSystemSolution()
        {
            // invertible
            var A1 = TriangularUpper.CreateFromArray(invertible);
            var b1 = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(invertible, vector10));
            var x1Expected = Vector.CreateFromArray(vector10);
            Vector x1Computed = A1.SolveLinearSystem(b1);
            comparer.AssertEqual(x1Expected, x1Computed);

            // singular
            var A2 = TriangularUpper.CreateFromArray(singular);
            var b2 = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(singular, vector10));
            var x2Expected = Vector.CreateFromArray(vector10);
            Vector x2Computed = A2.SolveLinearSystem(b2);
            Assert.False(comparer.AreEqual(x2Expected, x2Computed));

            // invertible - solve transposed (forward substitution)
            Matrix A3 = Matrix.CreateFromArray(invertible).Invert().Transpose();
            Vector x3Expected = A3 * b1;
            Vector x3Computed = A1.SolveLinearSystem(b1, true);
            comparer.AssertEqual(x3Expected, x3Computed);
        }

        [Fact]
        public static void TestTransposition()
        {
            // invertible
            var A1 = TriangularUpper.CreateFromArray(invertible);
            var A1TransposeExpected = MatrixOperations.Transpose(invertible);
            TriangularLower A1TransposeComputed = A1.Transpose(false);
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // singular
            var A2 = TriangularUpper.CreateFromArray(singular);
            var A2TransposeExpected = MatrixOperations.Transpose(singular);
            TriangularLower A2TransposeComputed = A2.Transpose(false);
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
