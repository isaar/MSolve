using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing
{
    class TestSuiteSparse
    {
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_create_handle")]
        public static extern IntPtr CreateHandle();

        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_destroy_handle")]
        public static extern void DestroyHandle(ref IntPtr handle);

        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_destroy_factor")]
        public static extern void DestroyFactor(ref IntPtr factorizedMatrix, IntPtr handle);

        /// <summary>
        /// Factorize a symmetric matrix using cholesky algorithm. The matrix is in csc form, with only the upper triangle stored.
        /// Returns: a) -1 when the factorization is successfull, 
        /// b) i>0 when the matrix is not positive definite with i being the column at which failure occured(0 based indexing), 
        /// c) -2 for other failures, such as a too large matrix
        /// </summary>
        /// <param name="order"></param>
        /// <param name="nnz"></param>
        /// <param name="values"></param>
        /// <param name="rowIndices"></param>
        /// <param name="colOffsets"></param>
        /// <param name="factorizedMatrix"></param>
        /// <param name="handle"></param>
        /// <returns></returns>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_factorize_cscupper")]
        public static extern int FactorizeCSCUpper(int order, int nnz, double[] values, int[] rowIndices, int[] colOffsets, out IntPtr factorizedMatrix, IntPtr handle);

        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_solve")]
        public static extern void Solve(int order, IntPtr factorizedMatrix, double[] rhs, double[] outSolution, IntPtr handle);

        public static void Run()
        {
            // Define linear system
            const int n = 4;
            const int nnz = 7;
            int[] col_offsets = new int[n + 1] { 0, 1, 2, 5, nnz };
            int[] row_indices = new int[nnz] { 0, 1, 0, 1, 2, 1, 3 };
            double[] values = new double[nnz] { 4.0, 10.0, 2.0, 1.0, 8.0, 3.0, 9.0 };
            double[] rhs = new double[n] { 6.0, 14.0, 11.0, 12.0 };
            double[] solution = new double[n];

            IntPtr handle = CreateHandle();
            IntPtr factor;
            int status = FactorizeCSCUpper(n, nnz, values, row_indices, col_offsets, out factor, handle);
            if (status != -1)
            {
                Console.WriteLine("Factorization failed");
            }
            Solve(n, factor, rhs, solution, handle);
            DestroyFactor(ref factor, handle);
            DestroyHandle(ref handle);

            //Check solution
            double[] expected = new double[n] { 1.0, 1.0, 1.0, 1.0 };
            Comparer comparer = new Comparer(Comparer.PrintMode.Never, 1e-6);
            if (comparer.AreEqual(expected, solution)) Console.WriteLine("The linear system has been solved correctly.");
            else Console.WriteLine("ERROR in solving the linear system.");
            Console.Write("expected solution = ");
            VectorMKL.CreateFromArray(expected).WriteToConsole();
            Console.Write("computed solution = ");
            VectorMKL.CreateFromArray(solution).WriteToConsole();
        }
    }
}
