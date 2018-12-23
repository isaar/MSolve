using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Factorizations
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
        //private static void TestFactorization()
        //{
        //    TestSettings.RunMultiproviderTest(providers, delegate () {
        //    
        //    int n = SparsePosDef10by10.order;
        //    var skyline = SkylineMatrix.CreateFromArrays(n, SparsePosDef10by10.skylineValues,
        //        SparsePosDef10by10.skylineDiagOffsets, true, true);
        //    var expectedU = Matrix.CreateFromArray(SparsePosDef10by10.choleskyU);
        //    CholeskySkyline factorization = skyline.FactorCholesky(false);
        //    TriangularUpper computedU = factorization.GetFactorU(); 
        //    comparer.AssertEqual(expectedU, computedU);
        //    });
        //}

        //TODO: Skyline operations are not part of the MKL SparseBLAS provider yet. There are only provided by the managed provider 
        //[Theory]
        //[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        [Fact]
        private static void TestMultipleSystemsSolution(/*LinearAlgebraProviderChoice providers*/)
        {
            //TestSettings.RunMultiproviderTest(providers, delegate () {
            //

            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order, SparsePosDef10by10.SkylineValues,
                 SparsePosDef10by10.SkylineDiagOffsets, true, true);
            CholeskySkyline factor = skyline.FactorCholesky(false);
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
            CholeskySkyline factor = skyline.FactorCholesky(false);
            Vector xComputed = factor.SolveLinearSystem(b);
            comparer.AssertEqual(xExpected, xComputed);

            //});
        }
    }
}
