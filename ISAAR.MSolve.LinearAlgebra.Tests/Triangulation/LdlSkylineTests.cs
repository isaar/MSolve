using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Triangulation
{
    /// <summary>
    /// Tests for <see cref="LdlSkyline"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class LdlSkylineTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        //TODO: Skyline operations are not part of the MKL SparseBLAS provider yet. They are only provided by the managed provider 
        //[Theory]
        //[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        [Fact]
        private static void TestMultipleSystemsSolution(/*LinearAlgebraProviderChoice providers*/)
        {
            //TestSettings.RunMultiproviderTest(providers, delegate () {
            //

            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order, SparsePosDef10by10.SkylineValues,
                 SparsePosDef10by10.SkylineDiagOffsets, true, true);
            LdlSkyline factor = skyline.FactorLdl(false);
            var identity = Matrix.CreateIdentity(SparsePosDef10by10.Order);
            var inverse = Matrix.CreateZero(SparsePosDef10by10.Order, SparsePosDef10by10.Order);
            factor.SolveLinearSystems(identity, inverse);

            var matrixTimesInverse = MatrixOperations.MatrixTimesMatrix(SparsePosDef10by10.Matrix, inverse.CopyToArray2D());
            comparer.AssertEqual(identity.CopyToArray2D(), matrixTimesInverse);

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
            LdlSkyline factor = skyline.FactorLdl(false);
            Vector xComputed = factor.SolveLinearSystem(b);
            comparer.AssertEqual(xExpected, xComputed);

            //});
        }
    }
}
