using System;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
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

        [Fact]
        private static void TestCholeskySolver()
        {
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
            IntPtr handle = SuiteSparseUtilities.CreateCommon(0, 0);
            int status = SuiteSparseUtilities.FactorizeCSCUpper(n, nnz, values, rowIndices, colOffsets, out IntPtr factor, handle);
            Assert.True(status == -1);
            int nnzFactor = SuiteSparseUtilities.GetFactorNonZeros(factor);
            Console.WriteLine($"Before factorization: nnz = {nnz}");
            Console.WriteLine($"After factorization: nnz = {nnzFactor}");
            
            SuiteSparseUtilities.Solve(0, n, 1, factor, rhs, solution, handle);
            comparer.AssertEqual(solutionExpected, solution);

            SuiteSparseUtilities.DestroyFactor(ref factor, handle);
            SuiteSparseUtilities.DestroyCommon(ref handle);
        }

        [Fact]
        private static void TestReordering()
        {
            int order = SparseSymm5by5.order;
            int[] rowIndices = SparseSymm5by5.cscRowIndices;
            int[] colOffsets = SparseSymm5by5.cscColOffsets;
            int[] permutation = new int[order];
            IntPtr common = SuiteSparseUtilities.CreateCommon(0, 0);
            int status = SuiteSparseUtilities.ReorderAMDUpper(order, rowIndices.Length, rowIndices, colOffsets, permutation,
                out int factorNNZ, common);
            Assert.True(status == 1, "SuiteSparse reordering failed. A possible reason is the lack of enough available memory");
            comparer.AssertEqual(SparseSymm5by5.matlabPermutationAMD, permutation);
           
            SuiteSparseUtilities.DestroyCommon(ref common);
        }
    }
}
