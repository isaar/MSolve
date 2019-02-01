using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Triangulation
{
    /// <summary>
    /// Tests for <see cref="CholeskyCSparseNet"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CholeskyCSparseNetTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestSystemSolution()
        {
            int order = SparsePosDef10by10.Order;
            var skyline = SkylineMatrix.CreateFromArrays(order, SparsePosDef10by10.SkylineValues,
                SparsePosDef10by10.SkylineDiagOffsets, true, true);
            var dok = DokSymmetric.CreateFromSparseMatrix(skyline);
            var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

            (double[] cscValues, int[] cscRowIndices, int[] cscColOffsets) = dok.BuildSymmetricCscArrays(true);
            var factor = CholeskyCSparseNet.Factorize(order, cscValues.Length, cscValues, cscRowIndices, cscColOffsets);
            Vector xComputed = factor.SolveLinearSystem(b);
            comparer.AssertEqual(xExpected, xComputed);
        }
    }
}
