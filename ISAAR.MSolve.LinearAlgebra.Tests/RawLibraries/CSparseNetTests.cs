using System;
using System.IO;
using CSparse;
using CSparse.Double.Factorization;
using CSparse.IO;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

//TODO: add test for rectangular matrix
namespace ISAAR.MSolve.LinearAlgebra.Tests.RawLibraries
{
    public static class CSparseNetTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        private static string PosDefMatrixPath => Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparsePosDef40by40MatrixMarket.mat";

        [Fact]
        public static void TestLinearSystemSolution()
        {
            // Load matrix from a file.
            var A = MatrixMarketReader.ReadMatrix<double>(PosDefMatrixPath);

            int m = A.RowCount;
            int n = A.ColumnCount;
            Assert.True(m == n);

            // Create test data.
            var xExpected = CSparse.Double.Vector.Create(n, 1.0);
            var b = new double[m];

            // Compute right hand side vector b.
            A.Multiply(1.0, xExpected, 0.0, b);

            // Apply column ordering to A to reduce fill-in.
            var order = ColumnOrdering.MinimumDegreeAtPlusA;

            // Factorize
            var xComputed = new double[n];
            var chol = SparseCholesky.Create(A, order);

            // Solve Ax = b (do not overwrite x).
            chol.Solve(b, xComputed);

            // Check the solution
            comparer.AssertEqual(xExpected, xComputed);


            // Compute residual b - Ax (do not overwrite b).
            var r = new double[m];
            Array.Copy(b, r, m);
            A.Multiply(-1.0, xComputed, 1.0, r);
        }

        private static void TestNormalEquationsSolution()
        {
            // Load a regtangular matrix from a file.
            string filePath = null;
            var A = MatrixMarketReader.ReadMatrix<double>(filePath);

            int m = A.RowCount;
            int n = A.ColumnCount;

            // Create test data.
            var x = CSparse.Double.Vector.Create(n, 1.0);
            var b = new double[m];

            // Compute right hand side vector b.
            A.Multiply(1.0, x, 0.0, b);

            // Apply column ordering to A to reduce fill-in.
            var order = ColumnOrdering.MinimumDegreeAtPlusA;

            if (m > n)
            {
                var At = A.Transpose();
                var AtA = At.Multiply(A);

                var c = CSparse.Double.Vector.Create(n, 0.0);

                At.Multiply(b, c);

                var chol = SparseCholesky.Create(AtA, order);

                // Solve normal equation A'Ax = A'b (overwrite x).
                chol.Solve(c, x);
            }
            else
            {
                Assert.True(false);
            }

            // Compute residual b - Ax (overwrite b).
            A.Multiply(-1.0, x, 1.0, b);
        }
    }
}
