using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing
{
    public static class BLASlvl1Tests
    {
        /// <summary>
        /// BLAS level 1 to test if the installation works correctly.
        /// </summary>
        public static void Run()
        {
            AddVectors();
            MultiplyVectors();           
        }

        private static void AddVectors()
        {
            int n = 5;
            double[] a = { 1, 2, 3, 4, 5 };
            double[] b = { 10, 20, 30, 40, 50 };

            double[] cMKL = new double[5];
            Array.Copy(b, cMKL, n);
            CBlas.Daxpy(a.Length, 1.0, ref a[0], 1, ref cMKL[0], 1);
            double[] cExpected = Utilities.MatrixOperations.AddArrays(a, b);

            double tol = 1e-10;
            bool correct = true;
            for (int i = 0; i < n; ++i)
            {
                if (Math.Abs(cMKL[i] - cExpected[i]) > tol)
                {
                    correct = false;
                    break;
                }
            }
            if (correct)
            {
                Console.WriteLine("Vector addition is correct.");
            }
            else
            {
                Console.WriteLine("Vector addition is NOT correct.");
            }
        }

        public static void MultiplyVectors()
        {
            int n = 5;
            double[] a = { 1, 2, 3, 4, 5 };
            double[] b = { 10, 20, 30, 40, 50 };

            double dotMKL = CBlas.Ddot(n, ref a[0], 1, ref b[0], 1);
            double dotExpected = Utilities.MatrixOperations.DotMultiplyArrays(a, b);

            double tol = 1e-10;
            if (Math.Abs(dotMKL - dotExpected) <= tol)
            {
                Console.WriteLine("Vector dot multiplication is correct.");
            }
            else
            {
                Console.WriteLine("Vector dot multiplication is NOT correct.");
            }
        }
    }
}
