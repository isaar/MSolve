using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.MatrixBuilders
{
    /// <summary>
    /// Tests for <see cref="DokColMajor"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class DokColMajorTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        private static DokColMajor CreateDok(double[,] matrix)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            var dok = DokColMajor.CreateEmpty(m, n);
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (matrix[i, j] != 0.0) dok[i, j] = matrix[i, j];
                }
            }
            return dok;
        }

        [Fact]
        private static void TestBuildCSR()
        {
            int m = SparseRectangular10by5.numRows;
            int n = SparseRectangular10by5.numCols;
            var dense = Matrix.CreateFromArray(SparseRectangular10by5.matrix);
            DokColMajor dok = CreateDok(SparseRectangular10by5.matrix);

            // CSR with sorted col indices of each row
            CscMatrix cscSorted = dok.BuildCscMatrix(true);
            comparer.AssertEqual(dense, cscSorted);

            // CSR without sorting
            CscMatrix cscUnsorted = dok.BuildCscMatrix(false);
            comparer.AssertEqual(dense, cscUnsorted);
        }

        [Fact]
        private static void TestIndexer()
        {
            int m = SparseRectangular10by5.numRows;
            int n = SparseRectangular10by5.numCols;
            Matrix dense = Matrix.CreateFromArray(SparseRectangular10by5.matrix);
            DokColMajor dok = CreateDok(SparseRectangular10by5.matrix);
            comparer.AssertEqual(dense, dok);
        }
    }
}
