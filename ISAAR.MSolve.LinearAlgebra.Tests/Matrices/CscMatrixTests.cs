using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="CscMatrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CscMatrixTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestEquality()
        {
            var full = Matrix.CreateFromArray(SparseRectangular10by5.matrix);
            var csc = CscMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.cscValues, SparseRectangular10by5.cscRowIndices, SparseRectangular10by5.cscColOffsets,
                true);
            Assert.True(csc.Equals(full));
        }

        [Fact]
        private static void TestMatrixCopy()
        {
            var full = Matrix.CreateFromArray(SparseRectangular10by5.matrix);
            var csc = CscMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.cscValues, SparseRectangular10by5.cscRowIndices, SparseRectangular10by5.cscColOffsets,
                true);
            comparer.AssertEqual(full, csc.CopyToFullMatrix());
        }

        [Fact]
        private static void TestMatrixMatrixMultiplication()
        {
            var matrix5x5 = Matrix.CreateFromArray(SquareInvertible10by10.matrix).GetSubmatrix(0, 5, 0, 5); //TODO: add a 5x5 matrix and its products
            var matrix10x10 = Matrix.CreateFromArray(SquareInvertible10by10.matrix);
            var ATimesMatrix5x5 = Matrix.CreateFromArray(
                MatrixOperations.MatrixTimesMatrix(SparseRectangular10by5.matrix, matrix5x5.CopyToArray2D()));
            var ATimesTransposeMatrix5x5 = Matrix.CreateFromArray(
                MatrixOperations.MatrixTimesMatrix(SparseRectangular10by5.matrix, matrix5x5.Transpose().CopyToArray2D()));
            var transposeATimesMatrix10x10 = Matrix.CreateFromArray(MatrixOperations.MatrixTimesMatrix(
                MatrixOperations.Transpose(SparseRectangular10by5.matrix), matrix10x10.CopyToArray2D()));
            var transposeATimesTransposeMatrix10x10 = Matrix.CreateFromArray(MatrixOperations.MatrixTimesMatrix(
                MatrixOperations.Transpose(SparseRectangular10by5.matrix), matrix10x10.Transpose().CopyToArray2D()));
            var matrix10x10TimesA = Matrix.CreateFromArray(
                MatrixOperations.MatrixTimesMatrix(matrix10x10.CopyToArray2D(), SparseRectangular10by5.matrix));
            var transposeMatrix10x10TimesA = Matrix.CreateFromArray(
                MatrixOperations.MatrixTimesMatrix(matrix10x10.Transpose().CopyToArray2D(), SparseRectangular10by5.matrix));
            var matrix5x5TimesTransposeA = Matrix.CreateFromArray(
                MatrixOperations.MatrixTimesMatrix(matrix5x5.CopyToArray2D(),
                MatrixOperations.Transpose(SparseRectangular10by5.matrix)));
            var transposeMatrix5x5TimesTransposeA = Matrix.CreateFromArray(
                MatrixOperations.MatrixTimesMatrix(matrix5x5.Transpose().CopyToArray2D(),
                MatrixOperations.Transpose(SparseRectangular10by5.matrix)));

            var A = CscMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.cscValues, SparseRectangular10by5.cscRowIndices, SparseRectangular10by5.cscColOffsets,
                true);

            // MultiplyRight()
            comparer.AssertEqual(ATimesMatrix5x5, A.MultiplyRight(matrix5x5, false, false));
            comparer.AssertEqual(ATimesTransposeMatrix5x5, A.MultiplyRight(matrix5x5, false, true));
            comparer.AssertEqual(transposeATimesMatrix10x10, A.MultiplyRight(matrix10x10, true, false));
            comparer.AssertEqual(transposeATimesTransposeMatrix10x10, A.MultiplyRight(matrix10x10, true, true));

            // MultiplyLeft()
            comparer.AssertEqual(matrix10x10TimesA, A.MultiplyLeft(matrix10x10, false, false));
            comparer.AssertEqual(transposeMatrix10x10TimesA, A.MultiplyLeft(matrix10x10, false, true));
            comparer.AssertEqual(matrix5x5TimesTransposeA, A.MultiplyLeft(matrix5x5, true, false));
            comparer.AssertEqual(transposeMatrix5x5TimesTransposeA, A.MultiplyLeft(matrix5x5, true, true));
        }

        [Fact]
        private static void TestMatrixVectorMultiplication()
        {
            // MultiplyRight() - untransposed 
            var A = CscMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.cscValues, SparseRectangular10by5.cscRowIndices, SparseRectangular10by5.cscColOffsets,
                true);
            var x5 = Vector.CreateFromArray(SparseRectangular10by5.lhs5);
            var b10Expected = Vector.CreateFromArray(SparseRectangular10by5.rhs10);
            Vector b10Computed = A.MultiplyRight(x5, false);
            comparer.AssertEqual(b10Expected, b10Computed);

            // MultiplyRight() - transposed
            var x10 = Vector.CreateFromArray(SparseRectangular10by5.lhs10);
            var b5Expected = Vector.CreateFromArray(SparseRectangular10by5.rhs5);
            Vector b5Computed = A.MultiplyRight(x10, true);
            comparer.AssertEqual(b5Expected, b5Computed);
        }

        [Fact]
        private static void TestTransposition()
        {
            var matrix = CscMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.cscValues, SparseRectangular10by5.cscRowIndices, SparseRectangular10by5.cscColOffsets,
                true);
            var transposeExpected = Matrix.CreateFromArray(MatrixOperations.Transpose(SparseRectangular10by5.matrix));

            // TransposeToCSC()
            CscMatrix transposeCsc = matrix.TransposeToCSC();
            comparer.AssertEqual(transposeExpected, transposeCsc);

            // TransposeToCSR()
            CsrMatrix transposeCsr = matrix.TransposeToCSR(true);
            comparer.AssertEqual(transposeExpected, transposeCsr);
        }
    }
}
