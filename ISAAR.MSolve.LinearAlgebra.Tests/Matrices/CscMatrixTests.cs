using ISAAR.MSolve.LinearAlgebra.Commons;
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
        private static void TestClear()
        {
            var zero = Matrix.CreateZero(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols);
            var csc = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                true);
            csc.Clear();
            comparer.AssertEqual(zero, csc);
        }

        [Fact]
        private static void TestEquality()
        {
            var full = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
            var csc = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                true);
            Assert.True(csc.Equals(full));
        }

        [Fact]
        private static void TestGetColumn()
        {
            var matrix = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                true);
            for (int j = 0; j < SparseRectangular10by5.NumCols; ++j)
            {
                Vector colExpected = DenseStrategies.GetColumn(matrix, j);
                Vector colComputed = matrix.GetColumn(j);
                comparer.AssertEqual(colExpected, colComputed);
            }
        }

        [Fact]
        private static void TestGetRow()
        {
            var matrix = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                true);
            for (int i = 0; i < SparseRectangular10by5.NumRows; ++i)
            {
                Vector rowExpected = DenseStrategies.GetRow(matrix, i);
                Vector rowComputed = matrix.GetRow(i);
                comparer.AssertEqual(rowExpected, rowComputed);
            }
        }

        [Fact]
        private static void TestMatrixCopy()
        {
            var full = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
            var csc = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                true);
            comparer.AssertEqual(full, csc.CopyToFullMatrix());
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixMatrixMultiplication(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var matrix5x5 = Matrix.CreateFromArray(SquareInvertible10by10.Matrix).GetSubmatrix(0, 5, 0, 5); //TODO: add a 5x5 matrix and its products
                var matrix10x10 = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
                var ATimesMatrix5x5 = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(SparseRectangular10by5.Matrix, matrix5x5.CopyToArray2D()));
                var ATimesTransposeMatrix5x5 = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(SparseRectangular10by5.Matrix, matrix5x5.Transpose().CopyToArray2D()));
                var transposeATimesMatrix10x10 = Matrix.CreateFromArray(MatrixOperations.MatrixTimesMatrix(
                    MatrixOperations.Transpose(SparseRectangular10by5.Matrix), matrix10x10.CopyToArray2D()));
                var transposeATimesTransposeMatrix10x10 = Matrix.CreateFromArray(MatrixOperations.MatrixTimesMatrix(
                    MatrixOperations.Transpose(SparseRectangular10by5.Matrix), matrix10x10.Transpose().CopyToArray2D()));
                var matrix10x10TimesA = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(matrix10x10.CopyToArray2D(), SparseRectangular10by5.Matrix));
                var transposeMatrix10x10TimesA = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(matrix10x10.Transpose().CopyToArray2D(), SparseRectangular10by5.Matrix));
                var matrix5x5TimesTransposeA = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(matrix5x5.CopyToArray2D(),
                    MatrixOperations.Transpose(SparseRectangular10by5.Matrix)));
                var transposeMatrix5x5TimesTransposeA = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(matrix5x5.Transpose().CopyToArray2D(),
                    MatrixOperations.Transpose(SparseRectangular10by5.Matrix)));

                var A = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                    SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
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
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixVectorMultiplication(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // MultiplyRight() - untransposed 
                var A = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                    SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                    true);
                var x5 = Vector.CreateFromArray(SparseRectangular10by5.Lhs5);
                var b10Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs10);
                Vector b10Computed = A.Multiply(x5, false);
                comparer.AssertEqual(b10Expected, b10Computed);

                // MultiplyRight() - transposed
                var x10 = Vector.CreateFromArray(SparseRectangular10by5.Lhs10);
                var b5Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs5);
                Vector b5Computed = A.Multiply(x10, true);
                comparer.AssertEqual(b5Expected, b5Computed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixVectorMultiplicationIntoResult(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // The result vectors will first be set to some non zero values to make sure that the result overwrites 
                // them instead of being added to them.

                // MultiplyIntoResult() - untransposed 
                var A = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                    SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                    true);
                var x5 = Vector.CreateFromArray(SparseRectangular10by5.Lhs5);
                var b10Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs10);
                Vector b10Computed = Vector.CreateWithValue(SparseRectangular10by5.NumRows, 1.0);
                //Vector bComputed = Vector.CreateZero(SparseRectangular10by5.NumRows);
                A.MultiplyIntoResult(x5, b10Computed, false);
                comparer.AssertEqual(b10Expected, b10Computed);

                // MultiplyIntoResult() - transposed
                var x10 = Vector.CreateFromArray(SparseRectangular10by5.Lhs10);
                var b5Expected = Vector.CreateFromArray(SparseRectangular10by5.Rhs5);
                Vector b5Computed = Vector.CreateWithValue(SparseRectangular10by5.NumCols, 1.0);
                A.MultiplyIntoResult(x10, b5Computed, true);
                comparer.AssertEqual(b5Expected, b5Computed);
            });
        }

        [Fact] //TODO: If the explicit transposition becomes abstracted behind a provider, then this should also be a Theory
        private static void TestTransposition()
        {
            var matrix = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                true);
            var transposeExpected = Matrix.CreateFromArray(MatrixOperations.Transpose(SparseRectangular10by5.Matrix));

            // TransposeToCSC()
            CscMatrix transposeCsc = matrix.TransposeToCSC();
            comparer.AssertEqual(transposeExpected, transposeCsc);

            // TransposeToCSR()
            CsrMatrix transposeCsr = matrix.TransposeToCSR(true);
            comparer.AssertEqual(transposeExpected, transposeCsr);
        }
    }
}
