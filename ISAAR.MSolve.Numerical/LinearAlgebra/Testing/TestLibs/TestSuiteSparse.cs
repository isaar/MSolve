using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
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
            IntPtr handle = SuiteSparseUtilities.CreateCommon();
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
            var matrixDOK = new DOKSymmetricColMajor(4);
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
