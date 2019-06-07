using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="SkylineMatrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SkylineMatrixTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestArrayCopy()
        {
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            comparer.AssertEqual(SparsePosDef10by10.Matrix, skyline.CopyToArray2D());
        }

        [Fact]
        private static void TestClear()
        {
            var zero = Matrix.CreateZero(SparsePosDef10by10.Order, SparsePosDef10by10.Order);
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            skyline.Clear();
            comparer.AssertEqual(zero, skyline);
        }

        [Fact]
        private static void TestConversionFullToSkyline()
        {
            var matrices = new List<double[,]>
            {
                DiagonalIndefinite.BuildeIndefiniteMatrix(20), GlobalMatrixAssembly.GlobalMatrix, SparsePosDef10by10.Matrix,
                SparseSymm5by5.Matrix, SymmPosDef10by10.Matrix, SymmSingular10by10.Matrix
            };

            foreach (double[,] matrix in matrices)
            {
                var full = Matrix.CreateFromArray(matrix);
                var skylineFromArray = SkylineMatrix.CreateFromArray(matrix);
                var skylineFromFull = SkylineMatrix.CreateFromMatrix(full);
                comparer.AssertEqual(full, skylineFromArray);
                comparer.AssertEqual(full, skylineFromFull);
            }
        }

        [Fact]
        private static void TestEquality()
        {
            var full = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order, 
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            Assert.True(skyline.Equals(full));
        }

        [Fact]
        private static void TestGetColumn()
        {
            var matrix = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            for (int j = 0; j < SparsePosDef10by10.Order; ++j)
            {
                Vector colExpected = DenseStrategies.GetColumn(matrix, j);
                Vector colComputed = matrix.GetColumn(j);
                comparer.AssertEqual(colExpected, colComputed);
            }
        }

        [Fact]
        private static void TestGetRow()
        {
            var matrix = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            for (int i = 0; i < SparsePosDef10by10.Order; ++i)
            {
                Vector rowExpected = DenseStrategies.GetRow(matrix, i);
                Vector rowComputed = matrix.GetRow(i);
                comparer.AssertEqual(rowExpected, rowComputed);
            }
        }

        [Fact]
        private static void TestGetSubmatrix()
        {
            // These are useful for debugging
            //string outputPath = @"C:\Users\Serafeim\Desktop\output.txt";
            //var writer = new LinearAlgebra.Output.FullMatrixWriter();

            var matrix = Matrix.CreateFromArray(
                MultiDiagonalMatrices.CreateSymmetric(100, new int[] { 2, 4, 8, 16, 32, 64 }));
            var matrixSky = SkylineMatrix.CreateFromMatrix(matrix);

            var indices = new int[] { 0, 2, 4, 6, 12, 24, 32, 50, 64, 80 };
            var indicesPerm = new int[] { 32, 80, 64, 0, 12, 24, 6, 50, 4, 2 };
            int[] rowIndices = indicesPerm;
            var colIndices = new int[] { 90, 10, 20, 60, 40, 50, 0, 70, 80, 30 };

            Matrix subMatrixFull = matrixSky.GetSubmatrixSymmetricFull(indices);
            SkylineMatrix subMatrixSky = matrixSky.GetSubmatrixSymmetricSkyline(indices);
            //writer.WriteToFile(subMatrixSky, outputPath, true);
            CscMatrix subMatrixCsc = matrixSky.GetSubmatrixCsc(indices, indices);

            Matrix subMatrixPermFull = matrixSky.GetSubmatrixSymmetricFull(indicesPerm);
            SkylineMatrix subMatrixPermSky = matrixSky.GetSubmatrixSymmetricSkyline(indicesPerm);
            CscMatrix subMatrixPermCsc = matrixSky.GetSubmatrixCsc(indicesPerm, indicesPerm);

            CscMatrix subMatrixRectCsc = matrixSky.GetSubmatrixCsc(rowIndices, colIndices);

            Matrix subMatrixExpected = matrix.GetSubmatrix(indices, indices);
            //writer.WriteToFile(subMatrixExpected, outputPath, true);
            Assert.True(subMatrixExpected.Equals(subMatrixFull));
            Assert.True(subMatrixExpected.Equals(subMatrixSky));
            Assert.True(subMatrixExpected.Equals(subMatrixCsc));

            Matrix subMatrixPermExpected = matrix.GetSubmatrix(indicesPerm, indicesPerm);
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermFull));
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermSky));
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermCsc));

            Matrix subMatrixRectExpected = matrix.GetSubmatrix(rowIndices, colIndices);
            Assert.True(subMatrixRectExpected.Equals(subMatrixRectCsc));
        }

        //TODO: Skyline operations are not part of the MKL SparseBLAS provider yet. There are only provided by the managed provider 
        //[Theory]
        //[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        [Fact]
        private static void TestMatrixVectorMultiplication(/*LinearAlgebraProviderChoice providers*/)
        {
            //TestSettings.RunMultiproviderTest(providers, delegate () {
            //

            var A = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            var x = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
            var bExpected = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
            IVector bComputed = A.Multiply(x, false);
            comparer.AssertEqual(bExpected, bComputed);

            //});
        }

        //TODO: Skyline operations are not part of the MKL SparseBLAS provider yet. There are only provided by the managed provider 
        //[Theory]
        //[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        [Fact]
        private static void TestMatrixVectorMultiplicationIntoResult(/*LinearAlgebraProviderChoice providers*/)
        {
            //TestSettings.RunMultiproviderTest(providers, delegate () {
            //

            var A = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            var x = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
            var bExpected = Vector.CreateFromArray(SparsePosDef10by10.Rhs);

            // The result vectors will first be set to some non zero values to make sure that the result overwrites 
            // them instead of being added to them.
            Vector bComputed = Vector.CreateWithValue(A.NumRows, 1.0);
            A.MultiplyIntoResult(x, bComputed, false);
            comparer.AssertEqual(bExpected, bComputed);

            //});
        }
    }
}
