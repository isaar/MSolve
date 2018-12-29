using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Factorizations
{
    /// <summary>
    /// Tests for <see cref="CholeskySuiteSparse"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CholeskySuiteSparseTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-12);

        [SkippableFact]
        private static void TestRowAddition()
        {
            Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

            Matrix original = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            Vector rhs = Vector.CreateFromArray(SparsePosDef10by10.Rhs);

            // Start the matrix as diagonal
            var matrixExpected = Matrix.CreateIdentity(original.NumColumns);
            var dok = DokSymmetric.CreateIdentity(SparsePosDef10by10.Order);
            var factor = CholeskySuiteSparse.Factorize(dok.BuildSymmetricCscMatrix(true), false);

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Reference solution
                #region minors are not positive definite this way
                //matrixExpected.SetRow(i, newRowVector);
                //matrixExpected.SetColumn(i, newRowVector);
                #endregion
                matrixExpected.SetSubmatrix(0, 0, original.GetSubmatrix(0, i + 1, 0, i + 1)); //this way they are
                //Console.WriteLine($"\nOnly dofs [0, {i}]");
                //matrixWriter.WriteToConsole(matrixExpected);

                // Update matrix
                Vector newRowVector = matrixExpected.GetRow(i);
                factor.AddRow(i, SparseVector.CreateFromDense(newRowVector));

                // Solve new linear system
                Vector solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                Vector solutionComputed = factor.SolveLinearSystem(rhs);
                comparer.AssertEqual(solutionExpected, solutionComputed);
            }
        }

        // will probably not work since the matrix will not always be positive definite
        //[SkippableFact]
        private static void TestRowAdditionReverse()
        {
            Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

            Matrix original = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            Vector rhs = Vector.CreateFromArray(SparsePosDef10by10.Rhs);

            // Start the matrix as diagonal
            var matrixExpected = Matrix.CreateIdentity(original.NumColumns);
            var dok = DokSymmetric.CreateIdentity(SparsePosDef10by10.Order);
            var factor = CholeskySuiteSparse.Factorize(dok.BuildSymmetricCscMatrix(true), false);

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Update matrix
                Vector newRowVector = original.GetRow(i);
                matrixExpected.SetSubrow(i, newRowVector);
                matrixExpected.SetSubcolumn(i, newRowVector);
                //Console.WriteLine($"\nOnly dofs [0, {i}]");
                factor.AddRow(i, SparseVector.CreateFromDense(newRowVector));

                // Solve new linear system
                Vector solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                Vector solutionComputed = factor.SolveLinearSystem(rhs);
                comparer.AssertEqual(solutionExpected, solutionComputed);
            }
        }

        [SkippableFact]
        private static void TestRowDeletion()
        {
            Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

            Matrix original = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            Vector rhs = Vector.CreateFromArray(SparsePosDef10by10.Rhs);

            // Start the matrix from the original
            var matrixExpected = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            var dok = DokSymmetric.CreateEmpty(SparsePosDef10by10.Order);
            for (int j = 0; j < matrixExpected.NumColumns; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    if (matrixExpected[i, j] != 0) dok[i, j] = matrixExpected[i, j];
                }
            }
            var factor = CholeskySuiteSparse.Factorize(dok.BuildSymmetricCscMatrix(true), false);

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Update matrix
                Vector identityRow = Vector.CreateZero(matrixExpected.NumColumns);
                identityRow[i] = 1.0;
                matrixExpected.SetSubrow(i, identityRow);
                matrixExpected.SetSubcolumn(i, identityRow);
                //Console.WriteLine($"\nOnly dofs [{i + 1}, 10)");
                factor.DeleteRow(i);

                // Solve new linear system
                Vector solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                Vector solutionComputed = factor.SolveLinearSystem(rhs);
                comparer.AssertEqual(solutionExpected, solutionComputed);
            }
        }

        [SkippableFact]
        private static void TestSystemSolution1()
        {
            Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

            // Define linear system
            var rhs = Vector.CreateFromArray(new double[] { 6.0, 14.0, 11.0, 12.0 });
            var solutionExpected = Vector.CreateFromArray(new double[] { 1.0, 1.0, 1.0, 1.0 });
            var matrixDOK = DokSymmetric.CreateEmpty(4);
            matrixDOK[0, 0] = 4.0; matrixDOK[0, 2] = 2.0;
            matrixDOK[1, 1] = 10.0; matrixDOK[1, 2] = 1.0; matrixDOK[1, 3] = 3.0;
            matrixDOK[2, 2] = 8.0;
            matrixDOK[3, 3] = 9.0;
            SymmetricCscMatrix matrixCSC = matrixDOK.BuildSymmetricCscMatrix(true);

            //const int n = 4;
            //const int nnz = 7;
            //int[] colOffsets = new int[n + 1] { 0, 1, 2, 5, nnz };
            //int[] rowIndices = new int[nnz] { 0, 1, 0, 1, 2, 1, 3 };
            //double[] values = new double[nnz] { 4.0, 10.0, 2.0, 1.0, 8.0, 3.0, 9.0 };
            //SymmetricCSC matrixCSC = new SymmetricCSC(values, rowIndices, colOffsets, false);

            //Solve it using SuiteSparse
            using (var factor = CholeskySuiteSparse.Factorize(matrixCSC, true))
            {
                Vector solution = factor.SolveLinearSystem(rhs);
                comparer.AssertEqual(solutionExpected, solution);
            }
        }

        [SkippableFact]
        private static void CheckSystemSolution2()
        {
            Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

            int order = SparsePosDef10by10.Order;

            // Build the matrices and right hand sides
            var dense = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            //var skyline = SkylineMatrix.CreateFromArrays(order, SparsePositiveDefinite.skylineValues, 
            //    SparsePositiveDefinite.skylineDiagOffsets, false);
            //var dok = DOKSymmetricColMajor.CreateFromSparseMatrix(skyline);
            var dok = DokSymmetric.CreateEmpty(order);
            for (int j = 0; j < order; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    if (dense[i, j] != 0) dok[i, j] = dense[i, j];
                }
            }
            Vector b = Vector.CreateFromArray(SparsePosDef10by10.Rhs);
            Matrix B = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);

            // Solve using dense algebra
            CholeskyFull chol = dense.FactorCholesky();
            Matrix U = chol.GetFactorU();
            Matrix L = U.Transpose();
            Vector xSolveExpect = chol.SolveLinearSystem(b);
            Matrix XSolveExpect = dense.Invert() * B;
            Vector xBackExpect = U.Invert() * b;
            Matrix XBackExpect = U.Invert() * B;
            Vector xForwardExpect = L.Invert() * b;
            Matrix XForwardExpect = L.Invert() * B;

            // Solve using SuiteSparse
            var (values, rowIndices, colOffsets) = dok.BuildSymmetricCscArrays(true);
            CholeskySuiteSparse factor = CholeskySuiteSparse.Factorize(order, values.Length, values, rowIndices, colOffsets,
                true);
            Vector xSolveComput = factor.SolveLinearSystem(b);
            Matrix XSolveComput = factor.SolveLinearSystems(B);
            Vector xBackComput = factor.BackSubstitution(b);
            Matrix XBackComput = factor.BackSubstitutions(B);
            Vector xForwardComput = factor.ForwardSubstitution(b);
            Matrix XForwardComput = factor.ForwardSubstitutions(B);
            Vector xSolveComput2 = factor.BackSubstitution(factor.ForwardSubstitution(b));

            // Check results
            comparer.AssertEqual(xSolveExpect, xSolveComput);
            comparer.AssertEqual(XSolveExpect, XSolveComput);
            comparer.AssertEqual(xBackExpect, xBackComput);
            comparer.AssertEqual(XBackExpect, XBackComput);
            comparer.AssertEqual(xForwardExpect, xForwardComput);
            comparer.AssertEqual(XForwardExpect, XForwardComput);
        }
    }
}
