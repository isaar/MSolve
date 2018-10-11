using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="Matrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MatrixTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        internal static void TestAddition()
        {
            var A1 = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            var A2 = Matrix.CreateFromArray(SymmPosDef10by10.matrix);
            var expected = Matrix.CreateFromArray(
                MatrixOperations.LinearCombination(1.0, SquareSingular10by10.matrix, 1.0, SymmPosDef10by10.matrix));
            
            // operator+
            comparer.AssertEqual(expected, A1 + A2);
        }

        [Fact]
        private static void TestEquality()
        {
            // Equals(SkylineMatrix)
            var full1 = Matrix.CreateFromArray(SparsePosDef10by10.matrix);
            var skyline1 = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.order,
                SparsePosDef10by10.skylineValues, SparsePosDef10by10.skylineDiagOffsets, true, true);
            Assert.True(full1.Equals(skyline1));

            // Equals(CsrMatrix)
            var full2 = Matrix.CreateFromArray(SparseRectangular10by5.matrix);
            var csr2 = CsrMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.csrValues, SparseRectangular10by5.csrColIndices, SparseRectangular10by5.csrRowOffsets,
                true);
            Assert.True(full2.Equals(csr2));

            // Equals(CscMatrix)
            var full3 = Matrix.CreateFromArray(SparseRectangular10by5.matrix);
            var csc3 = CscMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.cscValues, SparseRectangular10by5.cscRowIndices, SparseRectangular10by5.cscColOffsets,
                true);
            Assert.True(full3.Equals(csc3));
        }

        [Fact]
        private static void TestLinearCombination()
        {
            var A1 = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            double scalar1 = 2.0;
            var A2 = Matrix.CreateFromArray(SymmPosDef10by10.matrix);
            double scalar2 = 3.5;
            var expected = Matrix.CreateFromArray(
                MatrixOperations.LinearCombination(scalar1, SquareSingular10by10.matrix, scalar2, SymmPosDef10by10.matrix));

            // LinearCombination()
            comparer.AssertEqual(expected, A1.LinearCombination(scalar1, A2, scalar2));

            // LinearCombinationIntoThis()
            Matrix temp = A1.Copy();
            temp.LinearCombinationIntoThis(scalar1, A2, scalar2);
            comparer.AssertEqual(expected, temp);
        }

        [Fact]
        private static void TestMatrixMatrixMultiplication()
        {
            var A1 = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            var A2 = Matrix.CreateFromArray(RectangularFullRank10by5.matrix);
            var expectedA1TimesA2 = Matrix.CreateFromArray(
                MatrixOperations.MatrixTimesMatrix(SquareSingular10by10.matrix, RectangularFullRank10by5.matrix));
            var expectedTransposeA2TimesA1 = Matrix.CreateFromArray(
                MatrixOperations.MatrixTimesMatrix(
                    MatrixOperations.Transpose(RectangularFullRank10by5.matrix), SquareSingular10by10.matrix));
        
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
        private static void TestMatrixVectorMultiplication()
        {
            // rectangular 10-by-5
            var A1 = Matrix.CreateFromArray(RectangularFullRank10by5.matrix);
            var x1 = Vector.CreateFromArray(RectangularFullRank10by5.lhs5);
            var b1Expected = Vector.CreateFromArray(RectangularFullRank10by5.rhs10);
            Vector b1Computed = A1.MultiplyRight(x1, false);
            comparer.AssertEqual(b1Expected, b1Computed);

            // rectangular 5-by-10
            double[,] fullRank5by10 = MatrixOperations.Transpose(RectangularFullRank10by5.matrix);
            var x2 = Vector.CreateFromArray(RectangularFullRank10by5.lhs10);
            var b2Expected = Vector.CreateFromArray(RectangularFullRank10by5.rhs5);
            Vector b2Computed = A1.MultiplyRight(x2, true);
            comparer.AssertEqual(b2Expected, b2Computed);

            // square invertible 10-by-10
            var A3 = Matrix.CreateFromArray(SquareInvertible10by10.matrix);
            var x3 = Vector.CreateFromArray(SquareInvertible10by10.lhs);
            var b3Expected = Vector.CreateFromArray(SquareInvertible10by10.rhs);
            Vector b3Computed = A3.MultiplyRight(x3, false);
            comparer.AssertEqual(b3Expected, b3Computed);

            // square singular 10-by-10 (rank = 8)
            var A4 = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            var x4 = Vector.CreateFromArray(SquareSingular10by10.lhs);
            var b4Expected = Vector.CreateFromArray(SquareSingular10by10.rhs);
            Vector b4Computed = A4.MultiplyRight(x4, false);
            comparer.AssertEqual(b4Expected, b4Computed);

            // square singular 10-by-10 (rank = 9)
            var A5 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.matrix);
            var x5 = Vector.CreateFromArray(SquareSingularSingleDeficiency10by10.lhs);
            var b5Expected = Vector.CreateFromArray(SquareSingularSingleDeficiency10by10.rhs);
            Vector b5Computed = A5.MultiplyRight(x5, false);
            comparer.AssertEqual(b5Expected, b5Computed);
        }

        [Fact]
        private static void TestScaling()
        {
            var matrix = Matrix.CreateFromArray(RectangularFullRank10by5.matrix);
            double scalar = 5.0;
            var expected = Matrix.CreateFromArray(MatrixOperations.Scale(scalar, RectangularFullRank10by5.matrix));

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
        private static void TestSubtraction()
        {
            var A1 = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            var A2 = Matrix.CreateFromArray(SymmPosDef10by10.matrix);
            var expected = Matrix.CreateFromArray(
                MatrixOperations.LinearCombination(1.0, SquareSingular10by10.matrix, -1.0, SymmPosDef10by10.matrix));

            // operator+
            comparer.AssertEqual(expected, A1 - A2);
        }

        [Fact]
        private static void TestTransposition()
        {
            // square
            var A1 = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            var A1TransposeExpected = MatrixOperations.Transpose(SquareSingular10by10.matrix);
            Matrix A1TransposeComputed = A1.Transpose();
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // rectangular
            var A2 = Matrix.CreateFromArray(RectangularFullRank10by5.matrix);
            var A2TransposeExpected = MatrixOperations.Transpose(RectangularFullRank10by5.matrix);
            Matrix A2TransposeComputed = A2.Transpose();
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
