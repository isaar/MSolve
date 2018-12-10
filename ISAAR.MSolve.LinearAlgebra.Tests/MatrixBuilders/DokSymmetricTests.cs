using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.MatrixBuilders
{
    /// <summary>
    /// Tests for <see cref="DokSymmetric"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class DokSymmetricTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        private static DokSymmetric CreateDok(double[,] symmMatrix)
        {
            int n = symmMatrix.GetLength(0);
            var dok = DokSymmetric.CreateEmpty(n);
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    if (symmMatrix[i, j] != 0.0) dok[i, j] = symmMatrix[i, j];
                }
            }
            return dok;
        }

        [Fact]
        private static void TestAddSubmatrix()
        {
            var k1 = Matrix.CreateFromArray(GlobalMatrixAssembly.SubMatrix1);
            var k2 = Matrix.CreateFromArray(GlobalMatrixAssembly.SubMatrix2);
            var k3 = Matrix.CreateFromArray(GlobalMatrixAssembly.SubMatrix3);
            var expectedK = Matrix.CreateFromArray(GlobalMatrixAssembly.GlobalMatrix);

            var computedK = DokSymmetric.CreateEmpty(GlobalMatrixAssembly.GlobalOrder);
            computedK.AddSubmatrixSymmetric(k1, GlobalMatrixAssembly.IndicesDictionary1);
            computedK.AddSubmatrixSymmetric(k2, GlobalMatrixAssembly.IndicesDictionary2);
            computedK.AddSubmatrixSymmetric(k3, GlobalMatrixAssembly.IndicesDictionary3);

            comparer.AssertEqual(expectedK, computedK);
        }

        [Fact]
        private static void TestGetColumn()
        {
            Matrix dense = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            DokSymmetric dok = CreateDok(SparsePosDef10by10.Matrix);

            for (int j = 0; j < SparsePosDef10by10.Order; ++j)
            {
                comparer.AssertEqual(dense.GetColumn(j), dok.GetColumn(j)); //TODO: have hardcoded columns to compare against
            }
        }

        [Fact]
        private static void TestIndexer()
        {
            Matrix dense = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            DokSymmetric dok = CreateDok(SparsePosDef10by10.Matrix);
            comparer.AssertEqual(dense, dok);
        }
    }
}
