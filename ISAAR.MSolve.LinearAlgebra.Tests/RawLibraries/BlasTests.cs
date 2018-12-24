using System;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Providers.MKL;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;


namespace ISAAR.MSolve.LinearAlgebra.Tests.RawLibraries
{
    /// <summary>
    /// Tests for Intel MKL library's BLAS functions.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MklBlasTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [SkippableFact]
        private static void RunComputeNetExample()
        {
            Skip.IfNot(TestSettings.TestMkl, TestSettings.MessageWhenSkippingMKL);

            //Variable typeA is not used and it triggers a warning
#pragma warning disable CS0219

            const int GENERAL_MATRIX = 0;
            const int UPPER_MATRIX = 1;
            const int LOWER_MATRIX = -1;

            int m = 3, n = 2, i, j;
            int lda = 3, ldb = 3, ldc = 3;
            int rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc;
            float alpha = 0.5f, beta = 2.0f;
            float[] a, b, c;
            CBLAS_LAYOUT layout = CBLAS_LAYOUT.CblasRowMajor;
            CBLAS_SIDE side = CBLAS_SIDE.CblasLeft;
            CBLAS_UPLO uplo = CBLAS_UPLO.CblasUpper;

            int ma, na, typeA;
            if (side == CBLAS_SIDE.CblasLeft)
            {
                rmaxa = m + 1;
                cmaxa = m;
                ma = m;
                na = m;
            }
            else
            {
                rmaxa = n + 1;
                cmaxa = n;
                ma = n;
                na = n;
            }
            rmaxb = m + 1;
            cmaxb = n;
            rmaxc = m + 1;
            cmaxc = n;
            a = new float[rmaxa * cmaxa];
            b = new float[rmaxb * cmaxb];
            c = new float[rmaxc * cmaxc];
            if (layout == CBLAS_LAYOUT.CblasRowMajor)
            {
                lda = cmaxa;
                ldb = cmaxb;
                ldc = cmaxc;
            }
            else
            {
                lda = rmaxa;
                ldb = rmaxb;
                ldc = rmaxc;
            }
            if (uplo == CBLAS_UPLO.CblasUpper) typeA = UPPER_MATRIX;
            else typeA = LOWER_MATRIX;
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    a[i + j * lda] = 1.0f;
                }
            }
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    c[i + j * ldc] = 1.0f;
                    b[i + j * ldb] = 2.0f;
                }
            }
            CBlas.Ssymm(layout, side, uplo, m, n, alpha, ref a[0], lda, ref b[0], ldb, beta, ref c[0], ldc);
        }

        [SkippableFact]
        private static void TestAxpy()
        {
            Skip.IfNot(TestSettings.TestMkl, TestSettings.MessageWhenSkippingMKL);

            int n = 5;
            double[] a = { 1, 2, 3, 4, 5 };
            double[] b = { 10, 20, 30, 40, 50 };
            double[] cExpected = MatrixOperations.LinearCombination(1.0, a, 1.0, b);

            double[] cComputed = new double[5];
            Array.Copy(b, cComputed, n);
            CBlas.Daxpy(a.Length, 1.0, ref a[0], 1, ref cComputed[0], 1);

            comparer.AssertEqual(cExpected, cComputed);
        }

        [SkippableFact]
        private static void TestDdot()
        {
            Skip.IfNot(TestSettings.TestMkl, TestSettings.MessageWhenSkippingMKL);

            int n = 5;
            double[] a = { 1, 2, 3, 4, 5 };
            double[] b = { 10, 20, 30, 40, 50 };
            double dotExpected = MatrixOperations.DotProduct(a, b);
            double dotComputed = CBlas.Ddot(n, ref a[0], 1, ref b[0], 1);
            comparer.AssertEqual(dotExpected, dotComputed);
        }

        [SkippableFact]
        private static void TestDgemv()
        {
            Skip.IfNot(TestSettings.TestMkl, TestSettings.MessageWhenSkippingMKL);

            CBLAS_LAYOUT layout = CBLAS_LAYOUT.CblasColMajor;
            CBLAS_TRANSPOSE transA = CBLAS_TRANSPOSE.CblasNoTrans;
            double alpha = 1.0;
            double beta = 0.0;
            int incX = 1;
            int incY = 1;
            int m = RectangularFullRank10by5.NumRows;
            int n = RectangularFullRank10by5.NumCols;
            int ldA = m;
            double[] A = Conversions.Array2DToFullColMajor(RectangularFullRank10by5.Matrix);
            double[] X = RectangularFullRank10by5.Lhs5;
            double[] Y = new double[m];
            LapackePInvokes.Dgemv(layout, transA, m, n, alpha, A, ldA, X, incX, beta, Y, incY);

            comparer.AssertEqual(RectangularFullRank10by5.Rhs10, Y);
        }
    }
}
