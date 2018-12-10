using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Vectors
{
    /// <summary>
    /// Tests for <see cref="Vector"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class VectorTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-10);

        [Fact]
        private static void TestAddition()
        {
            var v1 = Vector.CreateFromArray(TestVectors.Vector1);
            var v2 = Vector.CreateFromArray(TestVectors.Vector2);
            var expected = Vector.CreateFromArray(TestVectors.Sum);

            // operator+
            comparer.AssertEqual(expected, v1 + v2);

            // AddIntoThis()
            var temp = Vector.CreateFromVector(v1);
            temp.AddIntoThis(v2);
            comparer.AssertEqual(expected, temp);
        }

        [Fact]
        private static void TestAxpy()
        {
            var v1 = Vector.CreateFromArray(TestVectors.Vector1);
            var v2 = Vector.CreateFromArray(TestVectors.Vector2);
            var expected = Vector.CreateFromArray(TestVectors.Vector1PlusVector2Times3);

            // Axpy()
            comparer.AssertEqual(expected, v1.Axpy(v2, TestVectors.Scalar2));

            // AxpyIntoThis
            var temp = Vector.CreateFromVector(v1);
            temp.AxpyIntoThis(v2, TestVectors.Scalar2);
            comparer.AssertEqual(expected, temp);
        }

        [Fact]
        private static void TestClear()
        {
            var zero = Vector.CreateZero(TestVectors.Vector1.Length);
            var vector = Vector.CreateFromArray(TestVectors.Vector1, true);
            vector.Clear();
            comparer.AssertEqual(zero, vector);
        }

        [Fact]
        private static void TestDotProduct()
        {
            var v1 = Vector.CreateFromArray(TestVectors.Vector1);
            var v2 = Vector.CreateFromArray(TestVectors.Vector2);

            // DotProduct()
            comparer.AssertEqual(TestVectors.DotProduct, v1.DotProduct(v2));

            // operator*
            comparer.AssertEqual(TestVectors.DotProduct, v1 * v2);
        }

        [Fact]
        private static void TestHadamardProduct()
        {
            var v1 = Vector.CreateFromArray(TestVectors.Vector1);
            var v2 = Vector.CreateFromArray(TestVectors.Vector2);
            var expected = Vector.CreateFromArray(TestVectors.HadamardProduct);

            // MultiplyPointwise()
            comparer.AssertEqual(expected, v1.MultiplyEntrywise(v2));

            // MultiplyPointwiseIntoThis()
            var temp = Vector.CreateFromVector(v1);
            temp.MultiplyEntrywiseIntoThis(v2);
            comparer.AssertEqual(expected, temp);
        }

        [Fact]
        private static void TestLinearCombination()
        {
            var v1 = Vector.CreateFromArray(TestVectors.Vector1);
            var v2 = Vector.CreateFromArray(TestVectors.Vector2);
            var expected = 2.5 * v1 + -3.5 * v2;
            var comparer = new MatrixComparer();

            // LinearCombination()
            comparer.AssertEqual(expected, v1.LinearCombination(2.5, v2, -3.5));

            // LinearCombinationIntoThis()
            var temp = Vector.CreateFromVector(v1);
            temp.LinearCombinationIntoThis(2.5, v2, -3.5);
            comparer.AssertEqual(expected, temp);
        }

        [Fact]
        private static void TestNorm2()
        {
            var vector = Vector.CreateFromArray(TestVectors.Vector1);

            // Norm2()
            comparer.AssertEqual(TestVectors.Norm2OfVector1, vector.Norm2());
        }

        [Fact]
        private static void TestScaling()
        {
            var vector = Vector.CreateFromArray(TestVectors.Vector1);
            var expected = Vector.CreateFromArray(TestVectors.Vector1Times2);

            // Scale()
            comparer.AssertEqual(expected, vector.Scale(2.0));

            // ScaleIntoThis()
            var temp = Vector.CreateFromVector(vector);
            temp.ScaleIntoThis(2.0);
            comparer.AssertEqual(expected, temp);

            // operator*
            comparer.AssertEqual(expected, 2.0 * vector);
        }

        [Fact]
        private static void TestSubtraction()
        {
            var v1 = Vector.CreateFromArray(TestVectors.Vector1);
            var v2 = Vector.CreateFromArray(TestVectors.Vector2);
            var expected = Vector.CreateFromArray(TestVectors.Difference);

            // operator-
            comparer.AssertEqual(expected, v1 - v2);

            // SubtractIntoThis()
            var temp = Vector.CreateFromVector(v1);
            temp.SubtractIntoThis(v2);
            comparer.AssertEqual(expected, temp);
        }
    }
}
