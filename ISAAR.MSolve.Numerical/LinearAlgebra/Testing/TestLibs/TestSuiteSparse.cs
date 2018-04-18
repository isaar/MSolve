using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestMatrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.SuiteSparse;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestLibs
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
            SuiteSparseUtilities.Solve(n, factor, rhs, solution, handle);
            SuiteSparseUtilities.DestroyFactor(ref factor, handle);
            SuiteSparseUtilities.DestroyCommon(ref handle);

            ProcessResult(solution);
        }

        public static void ExampleMatrixClasses()
        {
            // Define linear system
            var rhs = VectorMKL.CreateFromArray(new double[] { 6.0, 14.0, 11.0, 12.0 });
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
            using (CholeskySuiteSparse factor = matrixCSC.FactorCholesky())
            {
                VectorMKL solution = factor.SolveLinearSystem(rhs);
                ProcessResult(solution.CopyToArray());
            }
        }

        public static void CheckRowAddition()
        {
            FullMatrixWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 1, MaxIntegerDigits = 2 };
            FullVectorWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Array1DWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Comparer comparer = new Comparer(Comparer.PrintMode.Always);

            Matrix original = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            VectorMKL rhs = VectorMKL.CreateFromArray(SparsePositiveDefinite.rhs);
            Console.WriteLine("Full matrix: ");
            (new FullMatrixWriter(original)).WriteToConsole();

            // Start the matrix as diagonal
            var matrixExpected = Matrix.CreateIdentity(original.NumColumns);
            var dok = DOKSymmetricColMajor.CreateIdentity(SparsePositiveDefinite.order);
            CholeskySuiteSparse factor = dok.BuildSymmetricCSCMatrix(true).FactorCholesky();

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Update matrix
                matrixExpected.SetSubmatrix(0, 0, original.Slice(0, i + 1, 0, i + 1));
                Console.WriteLine($"\nOnly dofs [0, {i}]");
                (new FullMatrixWriter(matrixExpected)).WriteToConsole();
                VectorMKL newRowVector = matrixExpected.SliceRow(i);
                var newRow = new Dictionary<int, double>();
                for (int j = 0; j < newRowVector.Length; ++j)
                {
                    if (newRowVector[j] != 0) newRow.Add(j, newRowVector[j]);
                }
                factor.AddRow(i, SparseVector.CreateFromDictionary(matrixExpected.NumColumns, newRow));

                // Solve new linear system
                Console.WriteLine("\nCheck linear system solution");
                VectorMKL solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                VectorMKL solutionComputed = factor.SolveLinearSystem(rhs);
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
            VectorMKL rhs = VectorMKL.CreateFromArray(SparsePositiveDefinite.rhs);
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
            CholeskySuiteSparse factor = dok.BuildSymmetricCSCMatrix(true).FactorCholesky();

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Update matrix
                VectorMKL identityRow = VectorMKL.CreateZero(matrixExpected.NumColumns);
                identityRow[i] = 1.0;
                matrixExpected.SetRow(i, identityRow);
                matrixExpected.SetColumn(i, identityRow);
                Console.WriteLine($"\nOnly dofs [{i + 1}, 10)");
                (new FullMatrixWriter(matrixExpected)).WriteToConsole();
                factor.DeleteRow(i);

                // Solve new linear system
                Console.WriteLine("\nCheck linear system solution");
                VectorMKL solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                VectorMKL solutionComputed = factor.SolveLinearSystem(rhs);
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
