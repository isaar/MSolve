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
    /// Tests for <see cref="LUCSparseNet"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class LUCSparseNetTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestPosDefSystemSolution()
        {
            double pivotTolerance = 0.5;
            int order = SparsePosDef10by10.Order;
            var skyline = SkylineMatrix.CreateFromArrays(order, SparsePosDef10by10.SkylineValues,
                SparsePosDef10by10.SkylineDiagOffsets, true, true);
            var dok = DokColMajor.CreateFromSparseMatrix(skyline);
            var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);

            (double[] cscValues, int[] cscRowIndices, int[] cscColOffsets) = dok.BuildCscArrays(true);
            var factor = LUCSparseNet.Factorize(order, cscValues.Length, cscValues, cscRowIndices, cscColOffsets, pivotTolerance);
            Vector xComputed = factor.SolveLinearSystem(b);
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestSymmSystemSolution()
        {
            double pivotTolerance = 0.5;
            int order = SparseSymm5by5.Order;
            var csc = CscMatrix.CreateFromArrays(order, order, SparseSymm5by5.CscValues, SparseSymm5by5.CscRowIndices, 
                SparseSymm5by5.CscColOffsets, true);
            var xExpected = Vector.CreateWithValue(order, 1.0);
            var b = csc.Multiply(xExpected);

            var factor = LUCSparseNet.Factorize(order, csc.NumNonZeros, csc.RawValues, csc.RawRowIndices, csc.RawColOffsets,
                pivotTolerance);
            Vector xComputed = factor.SolveLinearSystem(b);
            comparer.AssertEqual(xExpected, xComputed);
        }

        [Fact]
        private static void TestNonsymmSystemSolution()
        {
            double pivotTolerance = 0.5;
            int order = SquareInvertible10by10.Order;
            var b = Vector.CreateFromArray(SquareInvertible10by10.Rhs);
            var xExpected = Vector.CreateFromArray(SquareInvertible10by10.Lhs);

            var dok = DokColMajor.CreateEmpty(order, order);
            for (int j = 0; j < order; ++j)
            {
                for (int i = 0; i < order; ++i) dok[i, j] = SquareInvertible10by10.Matrix[i, j];
            }

            (double[] cscValues, int[] cscRowIndices, int[] cscColOffsets) = dok.BuildCscArrays(true);
            var factor = LUCSparseNet.Factorize(order, cscValues.Length, cscValues, cscRowIndices, cscColOffsets, 
                pivotTolerance);
            Vector xComputed = factor.SolveLinearSystem(b);
            comparer.AssertEqual(xExpected, xComputed);
        }
    }
}
