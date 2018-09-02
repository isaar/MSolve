using ISAAR.MSolve.LinearAlgebra.Matrices;
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
    }
}
