using System;
using ISAAR.MSolve.LinearAlgebra.Providers.PInvoke;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.RawLibraries
{
    /// <summary>
    /// Tests for SuiteSparse functions.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SuiteSparseTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [SkippableFact]
        private static void TestCholeskySolver()
        {
            Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

            // Define linear system
            const int n = 4;
            const int nnz = 7;
            int[] colOffsets = new int[n + 1] { 0, 1, 2, 5, nnz };
            int[] rowIndices = new int[nnz] { 0, 1, 0, 1, 2, 1, 3 };
            double[] values = new double[nnz] { 4.0, 10.0, 2.0, 1.0, 8.0, 3.0, 9.0 };
            double[] rhs = new double[n] { 6.0, 14.0, 11.0, 12.0 };
            double[] solutionExpected = { 1.0, 1.0, 1.0, 1.0 };
            double[] solution = new double[n];

            // Solve it using SuiteSparse
            IntPtr handle = SuiteSparsePInvokes.CreateCommon(0, 0);
            int status = SuiteSparsePInvokes.FactorizeCSCUpper(n, nnz, values, rowIndices, colOffsets, out IntPtr factor, handle);
            Assert.True(status == -1);
            int nnzFactor = SuiteSparsePInvokes.GetFactorNonZeros(factor);
            Console.WriteLine($"Before factorization: nnz = {nnz}");
            Console.WriteLine($"After factorization: nnz = {nnzFactor}");
            
            SuiteSparsePInvokes.Solve(0, n, 1, factor, rhs, solution, handle);
            comparer.AssertEqual(solutionExpected, solution);

            SuiteSparsePInvokes.DestroyFactor(ref factor, handle);
            SuiteSparsePInvokes.DestroyCommon(ref handle);
        }

        [SkippableFact]
        private static void TestReordering()
        {
            Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

            int order = SparseSymm5by5.Order;
            int[] rowIndices = SparseSymm5by5.CscRowIndices;
            int[] colOffsets = SparseSymm5by5.CscColOffsets;
            int[] permutation = new int[order];
            IntPtr common = SuiteSparsePInvokes.CreateCommon(0, 0);
            int status = SuiteSparsePInvokes.ReorderAMDUpper(order, rowIndices.Length, rowIndices, colOffsets, permutation,
                out int factorNNZ, common);
            Assert.True(status == 1, "SuiteSparse reordering failed. A possible reason is the lack of enough available memory");
            comparer.AssertEqual(SparseSymm5by5.MatlabPermutationAMD, permutation);
           
            SuiteSparsePInvokes.DestroyCommon(ref common);
        }
    }
}
