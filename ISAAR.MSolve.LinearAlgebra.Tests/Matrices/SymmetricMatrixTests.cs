using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="SymmetricMatrix"/>
    /// </summary>
    public static class SymmetricMatrixTests
    {
        /// <summary>
        /// A symmetric 10-by-10 matrix that is positive definite.
        /// </summary>
        public static readonly double[,] symmPosDef10by10 = new double[,] {
            { 8.9156,    0.4590,    0.0588,    0.5776,    0.7118,    0.7423,    0.4389,    0.4353,    0.4929,    0.7223 },
            { 0.4590,    7.5366,    0.4276,    0.3282,    0.5277,    0.4274,    0.2498,    0.7622,    0.4987,    0.8953 },
            { 0.0588,    0.4276,    6.4145,    0.4144,    0.5954,    0.6196,    0.3257,    0.5084,    0.6342,    0.5270 },
            { 0.5776,    0.3282,    0.4144,    7.7576,    0.6609,    0.9212,    0.8040,    0.2146,    0.5077,    0.3928 },
            { 0.7118,    0.5277,    0.5954,    0.6609,    9.0882,    0.5096,    0.3434,    0.6139,    0.7590,    0.5276 },
            { 0.7423,    0.4274,    0.6196,    0.9212,    0.5096,   10.5523,    0.4885,    0.7861,    0.5294,    0.7933 },
            { 0.4389,    0.2498,    0.3257,    0.8040,    0.3434,    0.4885,    7.5277,    0.6574,    0.4020,    0.2763 },
            { 0.4353,    0.7622,    0.5084,    0.2146,    0.6139,    0.7861,    0.6574,    7.5142,    0.3528,    0.6254 },
            { 0.4929,    0.4987,    0.6342,    0.5077,    0.7590,    0.5294,    0.4020,    0.3528,    7.8487,    0.3035 },
            { 0.7223,    0.8953,    0.5270,    0.3928,    0.5276,    0.7933,    0.2763,    0.6254,    0.3035,    9.6122 }};

        /// <summary>
        /// A symmetric 10-by-10 matrix that is singular.
        /// </summary>
        public static readonly double[,] symmSingular10by10 = new double[,] {
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

        /// <summary>
        /// A vector with length = 10.
        /// </summary>
        public static readonly double[] vector10 = 
            { 2.6621, 3.5825, 0.8965, 1.6827, 0.9386, 1.6096, 2.0193, 2.7428, 0.2437, 2.7637 };

        private static readonly Comparer comparer = new Comparer(1E-13);

        [Fact]
        public static void TestArrayCopy()
        {
            // positive definite
            var A1 = SymmetricMatrix.CreateFromArray(symmPosDef10by10);
            comparer.AssertEqual(symmPosDef10by10, A1.CopyToArray2D());

            // singular
            var A2 = SymmetricMatrix.CreateFromArray(symmSingular10by10);
            comparer.AssertEqual(symmSingular10by10, A2.CopyToArray2D());
        }

        [Fact]
        public static void TestMatrixVectorMultiplication()
        {
            // invertible
            var A1 = SymmetricMatrix.CreateFromArray(symmPosDef10by10);
            var x1 = Vector.CreateFromArray(vector10);
            var b1Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(symmPosDef10by10, vector10));
            Vector b1Computed = A1.MultiplyRight(x1);
            comparer.AssertEqual(b1Expected, b1Computed);

            // singular
            var A2 = SymmetricMatrix.CreateFromArray(symmSingular10by10);
            var x2 = Vector.CreateFromArray(vector10);
            var b2Expected = Vector.CreateFromArray(MatrixOperations.MatrixTimesVector(symmSingular10by10, vector10));
            Vector b2Computed = A2.MultiplyRight(x1);
            comparer.AssertEqual(b2Expected, b2Computed);
        }

        [Fact]
        public static void TestTransposition()
        {
            // invertible
            var A1 = SymmetricMatrix.CreateFromArray(symmPosDef10by10);
            var A1TransposeExpected = MatrixOperations.Transpose(symmPosDef10by10);
            SymmetricMatrix A1TransposeComputed = A1.Transpose(false);
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // singular
            var A2 = SymmetricMatrix.CreateFromArray(symmSingular10by10);
            var A2TransposeExpected = MatrixOperations.Transpose(symmSingular10by10);
            SymmetricMatrix A2TransposeComputed = A2.Transpose(false);
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
