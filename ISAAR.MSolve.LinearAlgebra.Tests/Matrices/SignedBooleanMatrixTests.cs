using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="SignedBooleanMatrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SignedBooleanMatrixTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        internal static SignedBooleanMatrix CreateMatrix(double[,] matrix)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            var booleanMatrix = new SignedBooleanMatrix(m, n);
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (matrix[i, j] == 1.0) booleanMatrix.AddEntry(i, j, true);
                    else if (matrix[i, j] == -1.0) booleanMatrix.AddEntry(i, j, false);
                }
            }
            return booleanMatrix;
        }

        [Fact]
        private static void TestConstruction()
        {
            // matrix 1
            var A1Expected = Matrix.CreateFromArray(SignedBoolean5by10.A1);
            SignedBooleanMatrix A1Computed = CreateMatrix(SignedBoolean5by10.A1);
            comparer.AssertEqual(A1Expected, A1Computed);

            // matrix 2
            var A2Expected = Matrix.CreateFromArray(SignedBoolean5by10.A2);
            SignedBooleanMatrix A2Computed = CreateMatrix(SignedBoolean5by10.A2);
            comparer.AssertEqual(A2Expected, A2Computed);
        }

        [Fact]
        private static void TestMatrixVectorMultiplication()
        {
            // matrix 1 - untransposed
            SignedBooleanMatrix A1 = CreateMatrix(SignedBoolean5by10.A1);
            Vector x10 = Vector.CreateFromArray(SignedBoolean5by10.x10, true);
            Vector A1TimesX10Expected = Vector.CreateFromArray(SignedBoolean5by10.A1TimesX10, true);
            Vector A1TimesX10Computed = A1.MultiplyRight(x10, false);
            comparer.AssertEqual(A1TimesX10Expected, A1TimesX10Computed);

            // matrix 1 - transposed
            Vector x5 = Vector.CreateFromArray(SignedBoolean5by10.x5, true);
            Vector transpA1TimesX5Expected = Vector.CreateFromArray(SignedBoolean5by10.transpA1TimesX5, true);
            Vector transpA1TimesX5Computed = A1.MultiplyRight(x5, true);
            comparer.AssertEqual(transpA1TimesX5Expected, transpA1TimesX5Computed);

            // matrix 2 - untransposed
            SignedBooleanMatrix A2 = CreateMatrix(SignedBoolean5by10.A2);
            Vector A2TimesX10Expected = Vector.CreateFromArray(SignedBoolean5by10.A2TimesX10, true);
            Vector A2TimesX10Computed = A2.MultiplyRight(x10, false);
            comparer.AssertEqual(A2TimesX10Expected, A2TimesX10Computed);

            // matrix 2 - transposed
            Vector transpA2TimesX5Expected = Vector.CreateFromArray(SignedBoolean5by10.transpA2TimesX5, true);
            Vector transpA2TimesX5Computed = A2.MultiplyRight(x5, true);
            comparer.AssertEqual(transpA2TimesX5Expected, transpA2TimesX5Computed);
        }
    }
}
