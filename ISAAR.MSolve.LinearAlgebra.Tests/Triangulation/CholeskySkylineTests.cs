using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Triangulation
{
    /// <summary>
    /// Tests for <see cref="CholeskySkyline"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CholeskySkylineTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        //[Theory]
        //[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        [Fact]
        private static void TestFactorization()
        {
            //TestSettings.RunMultiproviderTest(providers, delegate () {
            //
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order, SparsePosDef10by10.SkylineValues,
                SparsePosDef10by10.SkylineDiagOffsets, true, true);
            var expectedU = Matrix.CreateFromArray(SparsePosDef10by10.CholeskyU);
            CholeskySkyline factorization = skyline.FactorCholesky(false);
            TriangularUpper computedU = factorization.GetFactorU();
            comparer.AssertEqual(expectedU, computedU);
            //});
        }

        //TODO: Skyline operations are not part of the MKL SparseBLAS provider yet. There are only provided by the managed provider 
        //[Theory]
        //[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        [Fact]
        private static void TestSystemSolution(/*LinearAlgebraProviderChoice providers*/)
        {
            //TestSettings.RunMultiproviderTest(providers, delegate () {
            //

            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order, SparsePosDef10by10.SkylineValues,
                 SparsePosDef10by10.SkylineDiagOffsets, true, true);
            var b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.Lhs);
            CholeskySkyline factor = skyline.FactorCholesky(false);
            Vector xComputed = factor.SolveLinearSystem(b);
            comparer.AssertEqual(xExpected, xComputed);

            //});
        }
    }
}
