using System;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Reordering
{
    /// <summary>
    /// Tests for <see cref="SparsityPatternSymmetric"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SparsityPatternSymmetricTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestConnectIndices()
        {
            int n = GlobalMatrixAssembly.globalOrder;
            var dense = Matrix.CreateFromArray(GlobalMatrixAssembly.globalMatrix);
            var pattern = SparsityPatternSymmetric.CreateEmpty(n);
            pattern.ConnectIndices(GlobalMatrixAssembly.globalIndices1, true);
            pattern.ConnectIndices(GlobalMatrixAssembly.globalIndices2, true);
            pattern.ConnectIndices(GlobalMatrixAssembly.globalIndices3, true);

            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    bool denseHasZero = dense[i, j] == 0.0;
                    bool patternHasZero = !pattern.IsNonZero(i, j);
                    Assert.True(patternHasZero == denseHasZero);
                }
            }
        }

        [Fact]
        private static void TestReorderingAMD()
        {
            var pattern = SparsityPatternSymmetric.CreateFromDense(Matrix.CreateFromArray(SparsePosDef10by10.matrix));
            var orderingAlg = new OrderingAmd();
            (int[] permutation, ReorderingStatistics stats) = orderingAlg.FindPermutation(pattern);
            comparer.AssertEqual(SparsePosDef10by10.matlabPermutationAMD, permutation);
        }

        [Fact]
        private static void TestReorderingCAMD()
        {
            int n = SparsePosDef10by10.order;
            var pattern = SparsityPatternSymmetric.CreateFromDense(Matrix.CreateFromArray(SparsePosDef10by10.matrix));
            var orderingAlg = new OrderingCamd();
            (int[] permutation, ReorderingStatistics stats) = 
                orderingAlg.FindPermutation(pattern, SparsePosDef10by10.constraintsCAMD);

            var originalDiagonal = new double[n];
            var permutedDiagonal = new double[n];
            for (int i = 0; i < n; ++i) originalDiagonal[i] = SparsePosDef10by10.matrix[i, i];
            for (int i = 0; i < n; ++i) permutedDiagonal[i] = originalDiagonal[permutation[i]];

            var writer = new Array1DWriter();
            Console.Write("Permutation (new-to-old): ");
            writer.WriteToConsole(permutation);
            Console.Write("Original diagonal: ");
            writer.WriteToConsole(originalDiagonal);
            Console.Write("Permuted diagonal: ");
            writer.WriteToConsole(permutedDiagonal);

            comparer.AssertEqual(SparsePosDef10by10.permutationCAMD, permutation);
        }
    }
}
