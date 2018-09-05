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

        //[Fact]
        //private static void TestFactorization()
        //{
        //    int n = SparsePosDef10by10.order;
        //    var skyline = SkylineMatrix.CreateFromArrays(n, SparsePosDef10by10.skylineValues,
        //        SparsePosDef10by10.skylineDiagOffsets, true, true);
        //    var expectedU = Matrix.CreateFromArray(SparsePosDef10by10.choleskyU);
        //    CholeskySkyline factorization = skyline.FactorCholesky(false);
        //    TriangularUpper computedU = factorization.GetFactorU(); 
        //    comparer.AssertEqual(expectedU, computedU);
        //}

        [Fact]
        private static void TestSystemSolution()
        {
            var skyline = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.order, SparsePosDef10by10.skylineValues,
                 SparsePosDef10by10.skylineDiagOffsets, true, true);
            var b = Vector.CreateFromArray(SparsePosDef10by10.rhs);
            var xExpected = Vector.CreateFromArray(SparsePosDef10by10.lhs);
            CholeskySkyline factor = skyline.FactorCholesky(false);
            Vector xComputed = factor.SolveLinearSystem(b);
            comparer.AssertEqual(xExpected, xComputed);
            
        }
    }
}
