using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.MKL;
using ISAAR.MSolve.LinearAlgebra.Commons;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestLibs
{
    public static class TestMkl
    {
        internal static void TestDgemv()
        {
            CBLAS_LAYOUT layout = CBLAS_LAYOUT.CblasColMajor;
            CBLAS_TRANSPOSE transA = CBLAS_TRANSPOSE.CblasNoTrans;
            double alpha = 1.0;
            double beta = 0.0;
            int incX = 1;
            int incY = 1;
            int m = RectangularFullColRank.numRows;
            int n = RectangularFullColRank.numCols;
            int ldA = m;
            double[] A = Conversions.Array2DToFullColMajor(RectangularFullColRank.matrix);
            double[] X = RectangularFullColRank.lhs5;
            double[] Y = new double[m];
            LAPACKE.Dgemv(layout, transA, m, n, alpha, A, ldA, X, incX, beta, Y, incY);

            Comparer comparer = new Comparer();
            if (comparer.AreEqual(RectangularFullColRank.rhs5, Y)) Console.WriteLine("CBLAS.Dgemv() was successful");
            else Console.WriteLine("CBLAS.Dgemv() failed");
        }

        internal static void TestDgetrf_Dgetrs()
        {
            int layout = LAPACKE.LAPACK_COL_MAJOR;
            char transA = LAPACKE.LAPACK_NO_TRANSPOSE;
            int m = SquareInvertible.order;
            int n = m;
            int minDim = m < n ? m : n;
            int ldA = m;
            double[] A = Conversions.Array2DToFullColMajor(SquareInvertible.matrix); // will be overwritten with LU
            int[] iPiv = new int[minDim];
            int nRhs = 1;
            double[] B = new double[n];
            Array.Copy(SquareInvertible.rhs, B, n); // will be overwritten with solution
            int ldB = n;

            int infoFact = MklUtilities.DefaultInfo;
            infoFact = LAPACKE.Dgetrf(layout, m, n, A, ldA, iPiv);
            Console.Write("LAPACKE.Dgetrf() result: ");
            if (infoFact == MklUtilities.DefaultInfo)
            {
                // first check the default info value, since it lies in the other intervals.
                // info == default => the MKL call did not succeed. 
                // info > 0 should not be returned at all by MKL, but it is here for completion.
                Console.WriteLine("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else if (infoFact < 0)
            {
                Console.WriteLine($"The {-infoFact}th parameter has an illegal value."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else if (infoFact > 0)
            {
                Console.WriteLine("The factorization has been completed, but U is singular."
                    + $" The first zero pivot is U[{infoFact - 1}, {infoFact - 1}] = 0.");
            }
            else Console.WriteLine("LAPACKE.Dgetrf() was successful");

            int infoSolve = MklUtilities.DefaultInfo;
            infoSolve = LAPACKE.Dgetrs(layout, transA, n, nRhs, A, ldA, iPiv, B, ldB);
            Console.Write("LAPACKE.Dgetrs() result: ");
            if (infoSolve == MklUtilities.DefaultInfo)
            {
                // first check the default info value, since it lies in the other intervals.
                // info == default => the MKL call did not succeed. 
                // info > 0 should not be returned at all by MKL, but it is here for completion.
                Console.WriteLine("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else if (infoSolve < 0)
            {
                Console.WriteLine($"The {-infoSolve}th parameter has an illegal value."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else
            {
                Comparer comparer = new Comparer();
                if (comparer.AreEqual(SquareInvertible.lhs, B)) Console.WriteLine("LAPACKE.Dgetrs() was successful");
                else Console.WriteLine("LAPACKE.Dgetrs() failed");
            }

        }
    }
}
