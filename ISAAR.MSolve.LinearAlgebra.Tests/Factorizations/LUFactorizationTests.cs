using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Factorizations
{
    /// <summary>
    /// Tests for <see cref="LUFactorization"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class LUFactorizationTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Fact]
        private static void TestDeterminantInvertiblePositive()
        {
            // invertible (rank = 10) with positive det
            var A = Matrix.CreateFromArray(SquareInvertible10by10.matrix);
            LUFactorization factorization = A.FactorLU();
            double det = factorization.CalcDeterminant();
            comparer.AssertEqual(SquareInvertible10by10.determinant, det);
        }

        [Fact]
        private static void TestDeterminantInvertibleNegative()
        {
            // Switch 2 rows to make the det negative
            var A = Matrix.CreateFromArray(SquareInvertible10by10.matrix);
            Vector row0 = A.GetRow(0);
            Vector row9 = A.GetRow(9);
            A.SetSubrow(0, row9);
            A.SetSubrow(9, row0);
            LUFactorization factorization = A.FactorLU();
            double det = factorization.CalcDeterminant();
            comparer.AssertEqual(-SquareInvertible10by10.determinant, det);
        }

        [Fact]
        private static void TestDeterminantSingular1()
        {
            // singular (rank = 8)
            var A = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            LUFactorization factorization = A.FactorLU();
            double det = factorization.CalcDeterminant();
            comparer.AssertEqual(SquareSingular10by10.determinant, det);
        }

        [Fact]
        private static void TestDeterminantSingular2()
        {
            // singular (rank = 9)
            var A = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.matrix);
            LUFactorization factorization = A.FactorLU();
            double det = factorization.CalcDeterminant();
            comparer.AssertEqual(SquareSingularSingleDeficiency10by10.determinant, det);
        }

        [Fact]
        private static void TestFactorsLU()
        {
            // invertible (rank = 10)
            var A1 = Matrix.CreateFromArray(SquareInvertible10by10.matrix);
            var expectedL1 = Matrix.CreateFromArray(SquareInvertible10by10.factorL);
            var expectedU1 = Matrix.CreateFromArray(SquareInvertible10by10.factorU);
            LUFactorization factorization1 = A1.FactorLU();
            Matrix computedL1 = factorization1.GetFactorL();
            Matrix computedU1 = factorization1.GetFactorU();
            comparer.AssertEqual(expectedL1, computedL1);
            comparer.AssertEqual(expectedU1, computedU1);

            // singular (rank = 8)
            var A2 = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            var expectedL2 = Matrix.CreateFromArray(SquareSingular10by10.factorL);
            var expectedU2 = Matrix.CreateFromArray(SquareSingular10by10.factorU);
            LUFactorization factorization2 = A2.FactorLU();
            Matrix computedL2 = factorization2.GetFactorL();
            Matrix computedU2 = factorization2.GetFactorU();
            comparer.AssertEqual(expectedL2, computedL2);
            comparer.AssertEqual(expectedU2, computedU2);

            // singular (rank = 9)
            var A3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.matrix);
            var expectedL3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.factorL);
            var expectedU3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.factorU);
            LUFactorization factorization3 = A3.FactorLU();
            Matrix computedL3 = factorization3.GetFactorL();
            Matrix computedU3 = factorization3.GetFactorU();
            comparer.AssertEqual(expectedL3, computedL3);
            comparer.AssertEqual(expectedU3, computedU3);
        }

        [Fact]
        private static void TestInversion()
        {
            // invertible (rank = 10)
            var A1 = Matrix.CreateFromArray(SquareInvertible10by10.matrix);
            var inverseA1Expected = Matrix.CreateFromArray(SquareInvertible10by10.inverse);
            LUFactorization factorization1 = A1.FactorLU();
            Matrix inverseA1Computed = factorization1.Invert(true);
            comparer.AssertEqual(inverseA1Expected, inverseA1Computed);

            // singular (rank = 8)
            var A2 = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            LUFactorization factorization2 = A2.FactorLU();
            Assert.Throws<SingularMatrixException>(() => factorization2.Invert(true));

            // singular (rank = 9)
            var A3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.matrix);
            LUFactorization factorization3 = A3.FactorLU();
            Assert.Throws<SingularMatrixException>(() => factorization3.Invert(true));
        }

        [Fact]
        private static void TestSystemSolution()
        {
            // invertible (rank = 10)
            var A1 = Matrix.CreateFromArray(SquareInvertible10by10.matrix);
            var b1 = Vector.CreateFromArray(SquareInvertible10by10.rhs);
            var x1Expected = Vector.CreateFromArray(SquareInvertible10by10.lhs);
            LUFactorization factorization1 = A1.FactorLU();
            Vector x1Computed = factorization1.SolveLinearSystem(b1);
            comparer.AssertEqual(x1Expected, x1Computed);

            // singular (rank = 8)
            var A2 = Matrix.CreateFromArray(SquareSingular10by10.matrix);
            var b2 = Vector.CreateFromArray(SquareSingular10by10.rhs);
            LUFactorization factorization2 = A2.FactorLU();
            Assert.Throws<SingularMatrixException>(() => factorization2.SolveLinearSystem(b2));

            // singular (rank = 9)
            var A3 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.matrix);
            var b3 = Vector.CreateFromArray(SquareSingularSingleDeficiency10by10.rhs);
            LUFactorization factorization3 = A3.FactorLU();
            Assert.Throws<SingularMatrixException>(() => factorization3.SolveLinearSystem(b3));
        }
    }
}
