using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Providers.MKL;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.RawLibraries
{
    /// <summary>
    /// Tests for Intel MKL library's LAPACK functions.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MklLapackTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [SkippableFact]
        internal static void TestDgetrf_Dgetrs()
        {
            Skip.IfNot(TestSettings.TestMkl, TestSettings.MessageWhenSkippingMKL);

            int layout = LapackePInvokes.LAPACK_COL_MAJOR;
            char transA = LapackePInvokes.LAPACK_NO_TRANSPOSE;
            int m = SquareInvertible10by10.Order;
            int n = m;
            int minDim = m < n ? m : n;
            int ldA = m;
            double[] A = Conversions.Array2DToFullColMajor(SquareInvertible10by10.Matrix); // will be overwritten with LU
            int[] iPiv = new int[minDim];
            int nRhs = 1;
            double[] B = new double[n];
            Array.Copy(SquareInvertible10by10.Rhs, B, n); // will be overwritten with solution
            int ldB = n;

            int infoFact = LapackUtilities.DefaultInfo;
            infoFact = LapackePInvokes.Dgetrf(layout, m, n, A, ldA, iPiv);
            Assert.True(infoFact == 0);
            //Console.Write("LAPACKE.Dgetrf() result: ");
            //if (infoFact == MklUtilities.DefaultInfo)
            //{
            //    // first check the default info value, since it lies in the other intervals.
            //    // info == default => the MKL call did not succeed. 
            //    // info > 0 should not be returned at all by MKL, but it is here for completion.
            //    Console.WriteLine("Something went wrong with the MKL call."
            //        + " Please contact the developer responsible for the linear algebra project.");
            //}
            //else if (infoFact < 0)
            //{
            //    Console.WriteLine($"The {-infoFact}th parameter has an illegal value."
            //        + " Please contact the developer responsible for the linear algebra project.");
            //}
            //else if (infoFact > 0)
            //{
            //    Console.WriteLine("The factorization has been completed, but U is singular."
            //        + $" The first zero pivot is U[{infoFact - 1}, {infoFact - 1}] = 0.");
            //}
            //else Console.WriteLine("LAPACKE.Dgetrf() was successful");

            int infoSolve = LapackUtilities.DefaultInfo;
            infoSolve = LapackePInvokes.Dgetrs(layout, transA, n, nRhs, A, ldA, iPiv, B, ldB);
            Assert.True(infoSolve == 0);
            //Console.Write("LAPACKE.Dgetrs() result: ");
            //if (infoSolve == MklUtilities.DefaultInfo)
            //{
            //    // first check the default info value, since it lies in the other intervals.
            //    // info == default => the MKL call did not succeed. 
            //    // info > 0 should not be returned at all by MKL, but it is here for completion.
            //    Console.WriteLine("Something went wrong with the MKL call."
            //        + " Please contact the developer responsible for the linear algebra project.");
            //}
            //else if (infoSolve < 0)
            //{
            //    Console.WriteLine($"The {-infoSolve}th parameter has an illegal value."
            //        + " Please contact the developer responsible for the linear algebra project.");
            //}
            comparer.AssertEqual(SquareInvertible10by10.Lhs, B);
        }
    }
}
