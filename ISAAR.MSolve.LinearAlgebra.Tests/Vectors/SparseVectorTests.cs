using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Vectors
{
    /// <summary>
    /// Tests for <see cref="SparseVector"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SparseVectorTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-10);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestAxpy(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var sparse = SparseVector.CreateFromArrays(SparseVector10.Length, SparseVector10.NonZeroValues,
                    SparseVector10.NonZeroIndices, true, false);
                var dense = Vector.CreateFromArray(SparseVector10.OtherVector);
                var expected = Vector.CreateFromArray(SparseVector10.OtherVectorPlusThisVectorTimes3);

                // Axpy()
                comparer.AssertEqual(expected, dense.Axpy(sparse, 3.0));

                // AxpyIntoThis
                var temp = Vector.CreateFromVector(dense);
                temp.AxpyIntoThis(sparse, 3.0);
                comparer.AssertEqual(expected, temp);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestDotProduct(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var sparse = SparseVector.CreateFromArrays(SparseVector10.Length, SparseVector10.NonZeroValues,
                    SparseVector10.NonZeroIndices, true, false);
                var dense = Vector.CreateFromArray(SparseVector10.OtherVector);

                // DotProduct() with dense vector
                comparer.AssertEqual(SparseVector10.DotThisTimesOther, sparse.DotProduct(dense));

                // DotProduct() with itsedlf
                comparer.AssertEqual(SparseVector10.Norm2OfThis * SparseVector10.Norm2OfThis, sparse.DotProduct(sparse));
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestNorm2(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var vector = SparseVector.CreateFromArrays(SparseVector10.Length, SparseVector10.NonZeroValues,
                    SparseVector10.NonZeroIndices, true, false);

                // Norm2()
                comparer.AssertEqual(SparseVector10.Norm2OfThis, vector.Norm2());
            });
        }
    }
}
