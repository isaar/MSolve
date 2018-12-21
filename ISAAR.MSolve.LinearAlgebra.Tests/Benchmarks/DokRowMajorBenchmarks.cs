using System;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Benchmarks
{
    /// <summary>
    /// Benchmarks for <see cref="DokRowMajor"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class DokRowMajorBenchmarks
    {
        public static void CompareDokCsrMultiplication()
        {
            int numRows = 100000;
            int numCols = 10000;
            DokRowMajor dok = RandomUtilities.CreateRandomSparseMatrix(numRows, numCols, 0.15);
            Vector lhs = RandomUtilities.CreateRandomVector(numCols);

            var watch = new Stopwatch();
            watch.Start();
            Vector dokTimesLhs = dok.MultiplyRight(lhs);
            watch.Stop();
            long dokTime = watch.ElapsedMilliseconds;

            LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.Managed;
            watch.Restart();
            Vector csrUnsortedTimesLhs = dok.BuildCsrMatrix(false).Multiply(lhs, false);
            watch.Stop();
            long csrUnsortedTime = watch.ElapsedMilliseconds;

            watch.Restart();
            Vector csrSortedTimesLhs = dok.BuildCsrMatrix(true).Multiply(lhs, false);
            watch.Stop();
            long csrSortedTime = watch.ElapsedMilliseconds;

            LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
            watch.Restart();
            Vector csrUnsortedMklTimesLhs = dok.BuildCsrMatrix(false).Multiply(lhs, false);
            watch.Stop();
            long csrUnsortedMklTime = watch.ElapsedMilliseconds;

            watch.Restart();
            Vector csrSortedMklTimesLhs = dok.BuildCsrMatrix(true).Multiply(lhs, false);
            watch.Stop();
            long csrSortedMklTime = watch.ElapsedMilliseconds;

            Console.WriteLine("Checking correctness and performance of DOK, sorted and unsorted CSR * vector:");
            Console.WriteLine($"DOK * vector: time = {dokTime} ms");
            Console.WriteLine($"Unsorted CSR * vector (C#): time = {csrUnsortedTime} ms");
            Console.WriteLine($"Sorted CSR * vector (C#): time = {csrSortedTime} ms");
            Console.WriteLine($"Unsorted CSR * vector (MKL): time = {csrUnsortedMklTime} ms");
            Console.WriteLine($"Sorted CSR * vector (MKL): time = {csrSortedMklTime} ms");

            var comparer = new ValueComparer(1e-12);
            double errorUnsorted = (csrUnsortedTimesLhs - dokTimesLhs).Norm2() / dokTimesLhs.Norm2();
            double errorSorted = (csrSortedTimesLhs - dokTimesLhs).Norm2() / dokTimesLhs.Norm2();
            double errorUnsortedMkl = (csrUnsortedMklTimesLhs - dokTimesLhs).Norm2() / dokTimesLhs.Norm2();
            double errorSortedMkl = (csrSortedMklTimesLhs - dokTimesLhs).Norm2() / dokTimesLhs.Norm2();
            Console.WriteLine("Multiplication DOK - unsorted CSR (C#): normalized error = " + errorUnsorted);
            Console.WriteLine("Multiplication DOK - sorted CSR (C#): normalized error = " + errorSorted);
            Console.WriteLine("Multiplication DOK - unsorted CSR (MKL): normalized error = " + errorUnsortedMkl);
            Console.WriteLine("Multiplication DOK - sorted CSR (MKL): normalized error = " + errorSortedMkl);
        }
    }
}
