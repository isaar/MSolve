using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="Matrix"/>
    /// </summary>
    public static class MatrixTests
    {
        /// <summary>
        /// A 10-by-5 rectangular matrix with full column rank (rank = number of columns). Its transpose will have full row rank.
        /// </summary>
        public static readonly double[,] fullRank10by5 = new double[,] {
            { 1.687338, 1.910625, 0.078085, 2.332444, 3.131424 },
            { 3.215230, 1.767679, 1.883659, 0.809185, 2.026953 },
            { 0.981285, 1.325970, 2.639193, 0.923645, 3.475809 },
            { 0.751857, 2.296387, 2.101308, 1.871156, 3.148936 },
            { 3.446088, 3.932793, 0.828995, 0.343996, 1.349877 },
            { 0.779645, 1.764859, 1.694237, 2.854230, 1.260338 },
            { 0.968975, 3.746153, 2.681011, 2.329230, 2.481479 },
            { 3.581092, 1.647008, 3.430712, 2.409152, 3.780291 },
            { 1.536053, 0.897209, 3.872541, 1.485733, 0.601071 },
            { 3.186167, 0.066731, 1.536828, 2.245892, 2.542407 }};

        /// <summary>
        /// A square 10-by-10 matrix that is invertible, but LU needs pivoting. 
        /// </summary>
        public static readonly double[,] invertible10by10 = new double[,] {
            { 1.0000,    2.6667,    7.0000,    5.0000,    2.5000,    9.0000,    6.0000,    2.2500,    4.0000,    3.0000 },
            { 0.0000,    0.3333,    1.0000,   -2.5000,    5.5000,   -7.7500,    2.0000,    5.7500,    3.0000,    2.0000 },
            { 2.0000,    2.0000,    4.0000,    9.0000,    1.7500,    5.0000,    3.5000,    2.0000,    6.0000,    8.0000 },
            { 1.0000,    1.0000,    3.0000,    5.0000,    1.7500,    6.0000,    4.0000,    7.0000,    6.0000,    9.0000 },
            { 5.0000,    0.6667,    4.5000,    5.0000,    0.2500,    2.3333,    1.5000,    1.2500,    9.0000,    0.7500 },
            { 0.3333,    2.0000,    4.0000,    5.0000,    5.0000,    9.0000,    2.5000,    4.0000,    1.0000,    3.0000 },
            { 0.2500,    6.0000,    8.0000,    0.7500,    1.0000,    3.0000,    8.0000,    1.2500,    1.0000,    2.0000 },
            { 1.0000,    1.0000,    1.0000,    7.0000,    1.0000,    5.0000,    2.0000,    7.0000,    1.0000,    2.0000 },
            { 3.5000,    9.0000,    9.0000,    2.0000,    8.0000,    4.0000,    2.5000,    6.0000,    7.0000,    9.0000 },
            { 3.0000,    7.0000,    9.0000,    4.0000,    2.6667,    9.0000,    8.0000,    5.0000,    7.0000,    3.0000 }};

        /// <summary>
        /// A square 10-by-10 matrix that is singular with rank = 8. 
        /// </summary>
        public static readonly double[,] singular10by10 = new double[,] {
            { 1.0000,    2.6667,    7.0000,    5.0000,    2.5000,    9.0000,    6.0000,    2.2500,    4.0000,    3.0000 },
            { 0.0000,    0.3333,    1.0000,   -2.5000,    5.5000,   -7.7500,    2.0000,    5.7500,    3.0000,    2.0000 },
            { 2.0000,    2.0000,    4.0000,    9.0000,    1.7500,    5.0000,    3.5000,    2.0000,    6.0000,    8.0000 },
            { 1.0000,    1.0000,    3.0000,    5.0000,    1.7500,    6.0000,    4.0000,    7.0000,    6.0000,    9.0000 },
            { 5.0000,    0.6667,    4.5000,    5.0000,    0.2500,    2.3333,    1.5000,    1.2500,    9.0000,    0.7500 },
            { 0.3333,    2.0000,    4.0000,    5.0000,    5.0000,    9.0000,    2.5000,    4.0000,    1.0000,    3.0000 },
            { 0.2500,    6.0000,    8.0000,    0.7500,    1.0000,    3.0000,    8.0000,    1.2500,    1.0000,    2.0000 },
            { 1.0000,    1.0000,    1.0000,    7.0000,    1.0000,    5.0000,    2.0000,    7.0000,    1.0000,    2.0000 },
            { 1.0000,    1.0000,    1.0000,    7.0000,    1.0000,    5.0000,    2.0000,    7.0000,    1.0000,    2.0000 },
            { 1.0000,    1.0000,    1.0000,    7.0000,    1.0000,    5.0000,    2.0000,    7.0000,    1.0000,    2.0000 } };

        /// <summary>
        /// A square 10-by-10 matrix that is singular with rank = 9. Only the last pivot will be 0, which causes MKL to classify 
        /// it as invertible. 
        /// </summary>
        public static readonly double[,] singular10by10SingleDeficiency = new double[,] {
            { 5.0389,    0.2156,    9.4373,    8.2953,    4.0673,    3.8888,    4.5039,    2.7529,    5.7474,    1.1704 },
            { 6.4681,    5.5984,    5.4916,    8.4909,    6.6693,    4.5474,    2.0567,    7.1667,    3.2604,    8.1468 },
            { 3.0775,    3.0082,    7.2839,    3.7253,    9.3373,    2.4669,    8.9965,    2.8338,    4.5642,    3.2486 },
            { 1.3872,    9.3941,    5.7676,    5.9318,    8.1095,    7.8442,    7.6259,    8.9620,    7.1380,    2.4623 },
            { 6.6851,    5.9753,    3.7231,    6.5385,    9.8797,    1.4888,    1.2281,    8.3437,    4.3851,    3.9582 },
            { 3.6246,    2.8662,    4.4653,    9.3350,    7.5675,    9.1371,    2.8495,    3.9003,    7.2086,    3.7569 },
            { 7.8811,    8.0082,    6.4630,    6.6846,    4.1705,    5.5828,    6.7323,    4.9790,    0.1861,    5.4655 },
            { 7.8030,    8.9611,    5.2120,    2.0678,    9.7179,    5.9887,    6.6428,    6.9481,    6.7478,    5.6192 },
            { 6.6851,    5.9753,    3.7231,    6.5385,    9.8797,    1.4888,    1.2281,    8.3437,    4.3851,    3.9582 },
            { 1.3350,    8.8402,    9.3713,    0.7205,    8.6415,    8.9971,    4.0732,    6.0963,    4.3782,    3.9813 }};

        /// <summary>
        /// A vector with length = 5.
        /// </summary>
        public static readonly double[] vector5 = { 3.4484, 1.9563, 2.7385, 4.2828, 5.3064 };

        /// <summary>
        /// A vector with length = 10.
        /// </summary>
        public static readonly double[] vector10 = 
            { 3.4484, 1.9563, 2.7385, 4.2828, 5.3064, 4.3251, 0.1117, 4.0487, 2.6311, 2.6269 };

        private static readonly Comparer comparer = new Comparer(1E-13);

        [Fact]
        public static void TestAddition()
        {
            double[,] posDef10by10 = SymmetricMatrixTests.symmPosDef10by10;
            var A1 = Matrix.CreateFromArray(singular10by10);
            var A2 = Matrix.CreateFromArray(posDef10by10);
            var expected = Matrix.CreateFromArray(MatrixOperations.LinearCombination(1.0, singular10by10, 1.0, posDef10by10));
            
            // operator+
            comparer.AssertEqual(expected, A1 + A2);
        }

        [Fact]
        public static void TestLinearCombination()
        {
            double[,] posDef10by10 = SymmetricMatrixTests.symmPosDef10by10;
            var A1 = Matrix.CreateFromArray(singular10by10);
            double scalar1 = 2.0;
            var A2 = Matrix.CreateFromArray(posDef10by10);
            double scalar2 = 3.5;
            var expected = Matrix.CreateFromArray(
                MatrixOperations.LinearCombination(scalar1, singular10by10, scalar2, posDef10by10));

            // LinearCombination()
            comparer.AssertEqual(expected, A1.LinearCombination(scalar1, A2, scalar2));

            // LinearCombinationIntoThis()
            Matrix temp = A1.Copy();
            temp.LinearCombinationIntoThis(scalar1, A2, scalar2);
            comparer.AssertEqual(expected, temp);
        }

        [Fact]
        public static void TestMatrixMatrixMultiplication()
        {
            var A1 = Matrix.CreateFromArray(singular10by10);
            var A2 = Matrix.CreateFromArray(fullRank10by5);
            var expectedA1TimesA2 = Matrix.CreateFromArray(MatrixOperations.MatrixTimesMatrix(singular10by10, fullRank10by5));
            var expectedTransposeA2TimesA1 = Matrix.CreateFromArray(
                MatrixOperations.MatrixTimesMatrix(MatrixOperations.Transpose(fullRank10by5), singular10by10));
        
            // MultiplyRight() without transposition
            comparer.AssertEqual(expectedA1TimesA2, A1.MultiplyRight(A2, false, false));

            // operator*
            comparer.AssertEqual(expectedA1TimesA2, A1 * A2);

            // MultiplyRight() with transposition
            comparer.AssertEqual(expectedTransposeA2TimesA1, A2.MultiplyRight(A1, true, false));

            // MultiplyRight() with incorrect dimensions
            Assert.Throws<NonMatchingDimensionsException>(() => A2.MultiplyRight(A1, false, false));
        }

        [Fact]
        public static void TestMatrixVectorMultiplication()
        {
            // rectangular 10-by-5
            var A1 = Matrix.CreateFromArray(fullRank10by5);
            var x1 = Vector.CreateFromArray(vector5);
            var b1Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(fullRank10by5, vector5));
            Vector b1Computed = A1.MultiplyRight(x1, false);
            comparer.AssertEqual(b1Expected, b1Computed);

            // rectangular 5-by-10
            double[,] fullRank5by10 = MatrixOperations.Transpose(fullRank10by5);
            var x2 = Vector.CreateFromArray(vector10);
            var b2Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(fullRank5by10, vector10));
            Vector b2Computed = A1.MultiplyRight(x2, true);
            comparer.AssertEqual(b2Expected, b2Computed);

            // square invertible 10-by-10
            var A3 = Matrix.CreateFromArray(invertible10by10);
            var x3 = Vector.CreateFromArray(vector10);
            var b3Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(invertible10by10, vector10));
            Vector b3Computed = A3.MultiplyRight(x3, false);
            comparer.AssertEqual(b3Expected, b3Computed);

            // square singular 10-by-10 (rank = 8)
            var A4 = Matrix.CreateFromArray(singular10by10);
            var x4 = Vector.CreateFromArray(vector10);
            var b4Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(singular10by10, vector10));
            Vector b4Computed = A4.MultiplyRight(x4, false);
            comparer.AssertEqual(b4Expected, b4Computed);

            // square singular 10-by-10 (rank = 9)
            var A5 = Matrix.CreateFromArray(singular10by10SingleDeficiency);
            var x5 = Vector.CreateFromArray(vector10);
            var b5Expected = Vector.CreateFromArray(
                MatrixOperations.MatrixTimesVector(singular10by10SingleDeficiency, vector10));
            Vector b5Computed = A5.MultiplyRight(x5, false);
            comparer.AssertEqual(b5Expected, b5Computed);
        }

        [Fact]
        public static void TestScaling()
        {
            var matrix = Matrix.CreateFromArray(fullRank10by5);
            double scalar = 5.0;
            var expected = Matrix.CreateFromArray(MatrixOperations.Scale(scalar, fullRank10by5));

            // Scale()
            comparer.AssertEqual(expected, matrix.Scale(scalar));

            // ScaleIntoThis()
            Matrix temp = matrix.Copy();
            temp.ScaleIntoThis(scalar);
            comparer.AssertEqual(expected, temp);

            // operator*
            comparer.AssertEqual(expected, scalar * matrix);
        }

        [Fact]
        public static void TestSubtraction()
        {
            double[,] posDef10by10 = SymmetricMatrixTests.symmPosDef10by10;
            var A1 = Matrix.CreateFromArray(singular10by10);
            var A2 = Matrix.CreateFromArray(posDef10by10);
            var expected = Matrix.CreateFromArray(MatrixOperations.LinearCombination(1.0, singular10by10, -1.0, posDef10by10));

            // operator+
            comparer.AssertEqual(expected, A1 - A2);
        }

        [Fact]
        public static void TestTransposition()
        {
            // square
            var A1 = Matrix.CreateFromArray(singular10by10);
            var A1TransposeExpected = MatrixOperations.Transpose(singular10by10);
            Matrix A1TransposeComputed = A1.Transpose();
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // rectangular
            var A2 = Matrix.CreateFromArray(fullRank10by5);
            var A2TransposeExpected = MatrixOperations.Transpose(fullRank10by5);
            Matrix A2TransposeComputed = A2.Transpose();
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
