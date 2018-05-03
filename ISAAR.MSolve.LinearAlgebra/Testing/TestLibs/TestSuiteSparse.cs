using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestLibs
{
    class TestSuiteSparse
    {
        public static void ExampleRawArrays()
        {
            // Define linear system
            const int n = 4;
            const int nnz = 7;
            int[] colOffsets = new int[n + 1] { 0, 1, 2, 5, nnz };
            int[] rowIndices = new int[nnz] { 0, 1, 0, 1, 2, 1, 3 };
            double[] values = new double[nnz] { 4.0, 10.0, 2.0, 1.0, 8.0, 3.0, 9.0 };
            double[] rhs = new double[n] { 6.0, 14.0, 11.0, 12.0 };
            double[] solution = new double[n];

            // Solve it using SuiteSparse
            IntPtr handle = SuiteSparseUtilities.CreateCommon(0, 0);
            int status = SuiteSparseUtilities.FactorizeCSCUpper(n, nnz, values, rowIndices, colOffsets, out IntPtr factor, handle);
            if (status != -1)
            {
                Console.WriteLine("Factorization failed");
                return;
            }
            else
            {
                int nnzFactor = SuiteSparseUtilities.GetFactorNonZeros(factor);
                Console.WriteLine($"Before factorization: nnz = {nnz}");
                Console.WriteLine($"After factorization: nnz = {nnzFactor}");
            }
            SuiteSparseUtilities.Solve(n, factor, rhs, solution, handle);
            SuiteSparseUtilities.DestroyFactor(ref factor, handle);
            SuiteSparseUtilities.DestroyCommon(ref handle);

            ProcessResult(solution);
        }

        public static void ExampleMatrixClasses()
        {
            // Define linear system
            var rhs = Vector.CreateFromArray(new double[] { 6.0, 14.0, 11.0, 12.0 });
            var matrixDOK = DOKSymmetricColMajor.CreateEmpty(4);
            matrixDOK[0, 0] = 4.0; matrixDOK[0, 2] = 2.0;
            matrixDOK[1, 1] = 10.0; matrixDOK[1, 2] = 1.0; matrixDOK[1, 3] = 3.0;
            matrixDOK[2, 2] = 8.0;
            matrixDOK[3, 3] = 9.0;
            SymmetricCSC matrixCSC = matrixDOK.BuildSymmetricCSCMatrix(true);

            //const int n = 4;
            //const int nnz = 7;
            //int[] colOffsets = new int[n + 1] { 0, 1, 2, 5, nnz };
            //int[] rowIndices = new int[nnz] { 0, 1, 0, 1, 2, 1, 3 };
            //double[] values = new double[nnz] { 4.0, 10.0, 2.0, 1.0, 8.0, 3.0, 9.0 };
            //SymmetricCSC matrixCSC = new SymmetricCSC(values, rowIndices, colOffsets, false);

            Console.WriteLine("Matrix DOK: ");
            (new FullMatrixWriter(matrixDOK)).WriteToConsole();
            Console.WriteLine("Matrix CSC: ");
            (new FullMatrixWriter(matrixCSC)).WriteToConsole();

            //Solve it using SuiteSparse
            using (CholeskySuiteSparse factor = matrixCSC.FactorCholesky(SuiteSparseOrdering.Natural))
            {
                Vector solution = factor.SolveLinearSystem(rhs);
                ProcessResult(solution.CopyToArray());
            }
        }

        public static void CheckReordering1()
        {
            int order = ReorderMatrix.order;
            int[] rowIndices = ReorderMatrix.cscRowIndices;
            int[] colOffsets = ReorderMatrix.cscColOffsets;
            int[] permutation = new int[order];
            IntPtr common = SuiteSparseUtilities.CreateCommon(0, 0);
            int status = SuiteSparseUtilities.ReorderAMDUpper(order, rowIndices.Length, rowIndices, colOffsets, permutation, 
                out int factorNNZ, common);
            if (status == 0)
                Console.WriteLine("SuiteSparse reordering failed. A possible reason is the lack of enough available memory");
            else
            {
                Comparer comparer = new Comparer();
                bool success = comparer.AreEqual(ReorderMatrix.matlabPermutationAMD, permutation);
                if (success) Console.WriteLine("AMD reordering was successful. The result is as expected.");
                else Console.WriteLine("SuiteSparse reordering returned, but the result is not as expected.");
            }
            SuiteSparseUtilities.DestroyCommon(ref common);
        }

        public static void CheckRowAddition()
        {
            FullMatrixWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 1, MaxIntegerDigits = 2 };
            FullVectorWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Array1DWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Comparer comparer = new Comparer(Comparer.PrintMode.Always);

            Matrix original = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            Vector rhs = Vector.CreateFromArray(SparsePositiveDefinite.rhs);
            Console.WriteLine("Full matrix: ");
            (new FullMatrixWriter(original)).WriteToConsole();

            // Start the matrix as diagonal
            var matrixExpected = Matrix.CreateIdentity(original.NumColumns);
            var dok = DOKSymmetricColMajor.CreateIdentity(SparsePositiveDefinite.order);
            CholeskySuiteSparse factor = dok.BuildSymmetricCSCMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural);

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Update matrix
                Vector newRowVector = original.SliceRow(i);
                #region minors are not positive definite this way
                //matrixExpected.SetRow(i, newRowVector);
                //matrixExpected.SetColumn(i, newRowVector);
                #endregion
                matrixExpected.SetSubmatrix(0, 0, original.Slice(0, i + 1, 0, i + 1)); //this way they are
                Console.WriteLine($"\nOnly dofs [0, {i}]");
                (new FullMatrixWriter(matrixExpected)).WriteToConsole();
                factor.AddRow(i, SparseVector.CreateFromDense(newRowVector));

                // Solve new linear system
                Console.WriteLine("\nCheck linear system solution");
                Vector solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                Vector solutionComputed = factor.SolveLinearSystem(rhs);
                comparer.CheckVectorEquality(solutionExpected, solutionComputed);
            }
        }

        // will probably not work since the matrix will not always be positive definite
        public static void CheckRowAdditionReverse()
        {
            FullMatrixWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 1, MaxIntegerDigits = 2 };
            FullVectorWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Array1DWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Comparer comparer = new Comparer(Comparer.PrintMode.Always);

            Matrix original = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            Vector rhs = Vector.CreateFromArray(SparsePositiveDefinite.rhs);
            Console.WriteLine("Full matrix: ");
            (new FullMatrixWriter(original)).WriteToConsole();

            // Start the matrix as diagonal
            var matrixExpected = Matrix.CreateIdentity(original.NumColumns);
            var dok = DOKSymmetricColMajor.CreateIdentity(SparsePositiveDefinite.order);
            CholeskySuiteSparse factor = dok.BuildSymmetricCSCMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural);

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Update matrix
                Vector newRowVector = original.SliceRow(i);
                matrixExpected.SetRow(i, newRowVector);
                matrixExpected.SetColumn(i, newRowVector);
                Console.WriteLine($"\nOnly dofs [0, {i}]");
                (new FullMatrixWriter(matrixExpected)).WriteToConsole();
                factor.AddRow(i, SparseVector.CreateFromDense(newRowVector));

                // Solve new linear system
                Console.WriteLine("\nCheck linear system solution");
                Vector solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                Vector solutionComputed = factor.SolveLinearSystem(rhs);
                comparer.CheckVectorEquality(solutionExpected, solutionComputed);
            }
        }

        public static void CheckRowDeletion()
        {
            FullMatrixWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 1, MaxIntegerDigits = 2 };
            FullVectorWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Array1DWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Comparer comparer = new Comparer(Comparer.PrintMode.Always);

            Matrix original = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            Vector rhs = Vector.CreateFromArray(SparsePositiveDefinite.rhs);
            Console.WriteLine("Full matrix: ");
            (new FullMatrixWriter(original)).WriteToConsole();

            // Start the matrix from the original
            var matrixExpected = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            var dok = DOKSymmetricColMajor.CreateEmpty(SparsePositiveDefinite.order);
            for (int j = 0; j < matrixExpected.NumColumns; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    if (matrixExpected[i, j] != 0) dok[i, j] = matrixExpected[i, j];
                }
            }
            CholeskySuiteSparse factor = dok.BuildSymmetricCSCMatrix(true).FactorCholesky(SuiteSparseOrdering.Natural);

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Update matrix
                Vector identityRow = Vector.CreateZero(matrixExpected.NumColumns);
                identityRow[i] = 1.0;
                matrixExpected.SetRow(i, identityRow);
                matrixExpected.SetColumn(i, identityRow);
                Console.WriteLine($"\nOnly dofs [{i + 1}, 10)");
                (new FullMatrixWriter(matrixExpected)).WriteToConsole();
                factor.DeleteRow(i);

                // Solve new linear system
                Console.WriteLine("\nCheck linear system solution");
                Vector solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                Vector solutionComputed = factor.SolveLinearSystem(rhs);
                comparer.CheckVectorEquality(solutionExpected, solutionComputed);
            }
        }

        private static void ProcessResult(double[] solution)
        {
            double[] expected = { 1.0, 1.0, 1.0, 1.0 };
            Comparer comparer = new Comparer(Comparer.PrintMode.Never, 1e-6);
            if (comparer.AreEqual(expected, solution)) Console.WriteLine("The linear system has been solved correctly.");
            else Console.WriteLine("ERROR in solving the linear system.");
            Console.Write("expected solution = ");
            (new Array1DWriter(expected)).WriteToConsole();
            Console.Write("computed solution = ");
            (new Array1DWriter(solution)).WriteToConsole();
        }
    }
}
