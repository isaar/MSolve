using System;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    class SparseRandomMatrix
    {
        public static void CompareDokCsrMultiplication()
        {
            int numRows = 100000;
            int numCols = 10000;
            DOKRowMajor dok = CreateRandomMatrix(numRows, numCols, 0.15);
            Vector lhs = CreateRandomVector(numCols);

            var watch = new Stopwatch();
            watch.Start();
            Vector dokTimesLhs = dok.MultiplyRight(lhs);
            watch.Stop();
            long dokTime = watch.ElapsedMilliseconds;

            CSRMatrix.UseMKL = false;
            watch.Restart();
            Vector csrUnsortedTimesLhs = dok.BuildCSRMatrix(false).MultiplyRight(lhs, false);
            watch.Stop();
            long csrUnsortedTime = watch.ElapsedMilliseconds;

            watch.Restart();
            Vector csrSortedTimesLhs = dok.BuildCSRMatrix(true).MultiplyRight(lhs, false);
            watch.Stop();
            long csrSortedTime = watch.ElapsedMilliseconds;

            CSRMatrix.UseMKL = true;
            watch.Restart();
            Vector csrUnsortedMklTimesLhs = dok.BuildCSRMatrix(false).MultiplyRight(lhs, false);
            watch.Stop();
            long csrUnsortedMklTime = watch.ElapsedMilliseconds;

            watch.Restart();
            Vector csrSortedMklTimesLhs = dok.BuildCSRMatrix(true).MultiplyRight(lhs, false);
            watch.Stop();
            long csrSortedMklTime = watch.ElapsedMilliseconds;

            Console.WriteLine("Checking correctness and performance of DOK, sorted and unsorted CSR * vector:");
            Console.WriteLine($"DOK * vector: time = {dokTime} ms");
            Console.WriteLine($"Unsorted CSR * vector (C#): time = {csrUnsortedTime} ms");
            Console.WriteLine($"Sorted CSR * vector (C#): time = {csrSortedTime} ms");
            Console.WriteLine($"Unsorted CSR * vector (MKL): time = {csrUnsortedMklTime} ms");
            Console.WriteLine($"Sorted CSR * vector (MKL): time = {csrSortedMklTime} ms");

            var comparer = new ValueComparer(1e-12);
            double errorUnsorted = csrUnsortedTimesLhs.Subtract(dokTimesLhs).Norm2() / dokTimesLhs.Norm2();
            double errorSorted = csrSortedTimesLhs.Subtract(dokTimesLhs).Norm2() / dokTimesLhs.Norm2();
            double errorUnsortedMkl = csrUnsortedMklTimesLhs.Subtract(dokTimesLhs).Norm2() / dokTimesLhs.Norm2();
            double errorSortedMkl = csrSortedMklTimesLhs.Subtract(dokTimesLhs).Norm2() / dokTimesLhs.Norm2();
            Console.WriteLine("Multiplication DOK - unsorted CSR (C#): normalized error = " + errorUnsorted);
            Console.WriteLine("Multiplication DOK - sorted CSR (C#): normalized error = " + errorSorted);
            Console.WriteLine("Multiplication DOK - unsorted CSR (MKL): normalized error = " + errorUnsortedMkl);
            Console.WriteLine("Multiplication DOK - sorted CSR (MKL): normalized error = " + errorSortedMkl);
        }

        public static DOKRowMajor CreateRandomMatrix(int numRows, int numCols, double nonZeroChance)
        {
            var rand = new Random();
            var dok = DOKRowMajor.CreateEmpty(numRows, numCols);
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numCols; ++j)
                {
                    if (rand.NextDouble() <= nonZeroChance)
                    {
                        dok[i, j] = rand.NextDouble();
                    }
                }
            }
            return dok;
        }

        public static Vector CreateRandomVector(int length)
        {
            var rand = new Random();
            var vector = new double[length];
            for (int i = 0; i < length; ++i)
            {
                vector[i] = rand.NextDouble();
            }
            return Vector.CreateFromArray(vector, false);
        }
    }
}
