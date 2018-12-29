using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Tests.GRACM;
using ISAAR.MSolve.XFEM.Tests.Khoei;
using ISAAR.MSolve.LinearAlgebra.Factorizations;

namespace ISAAR.MSolve.XFEM.Tests.DofOrdering
{
    class ReorderingTests
    {
        public static void Run()
        {
            Console.WriteLine("GRACM DCB: ");
            TestReordering(CreateDcbGRACM(), true, false);
            TestReordering(CreateDcbGRACM(), false, false);
            Console.WriteLine();

            Console.WriteLine("DCB 3x1: ");
            TestReordering(CreateDcb3x1(), true, false);
            TestReordering(CreateDcb3x1(), false, false);
            Console.WriteLine();

        }

        public static void TestAMDNonZeros(Model2D model, bool interleaved, bool printEnumerations)
        {
            TestReordering(model, interleaved, printEnumerations);
            TestReordering(model, interleaved, printEnumerations);
        }

        private static void TestReordering(Model2D model, bool interleaved, bool printEnumerations)
        {
            var watch = new Stopwatch();

            IDofOrderer unorderedDofs;
            string unorderedName;
            if (interleaved)
            {
                unorderedDofs = InterleavedDofOrderer.Create(model); ;
                unorderedName = "Interleaved";
            }
            else
            {
                unorderedDofs = SeparateDofOrderer.Create(model); ;
                unorderedName = "Separate";
            }

            // Before reorder
            var assembler = new GlobalDOKAssembler();
            (DokSymmetric Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model, unorderedDofs);
            Console.WriteLine($"CSC non zeros = {Kuu.CountNonZeros()}");

            watch.Start();
            using (var factor = CholeskySuiteSparse.Factorize(Kuu.BuildSymmetricCscMatrix(true), true))
            {
                watch.Stop();
                Console.WriteLine($"{unorderedName} ordering -> factorization: Non zeros = {factor.NumNonZerosUpper}"
                    + $" , time for factorization = {watch.ElapsedMilliseconds} ms");
            }

            // After reorder
            (int[] permutation, ReorderingStatistics stats) = (new OrderingAmdSuiteSparse()).FindPermutation(Kuu);
            Console.WriteLine($"{unorderedName} ordering -> AMD : Non zeros predicted by AMD = {stats.SupFactorizedNumNonZeros}");
            IDofOrderer reorderedDofs = unorderedDofs.DeepCopy();
            reorderedDofs.ReorderUnconstrainedDofs(permutation, false);
            (Kuu, Kuc) = assembler.BuildGlobalMatrix(model, reorderedDofs);

            watch.Restart();
            using (var factor = CholeskySuiteSparse.Factorize(Kuu.BuildSymmetricCscMatrix(true), true))
            {
                watch.Stop();
                Console.WriteLine($"{unorderedName} ordering -> AMD -> factorization: Non zeros = {factor.NumNonZerosUpper}"
                    + $" , time for factorization = {watch.ElapsedMilliseconds} ms");
            }

            // Deprecated: Let SuiteSparse handle the AMD ordering
            //(Kuu, Kuc) = assembler.BuildGlobalMatrix(model, unorderedDofs);
            //watch.Restart();
            //using (var factor = CholeskySuiteSparse.Factorize(Kuu.BuildSymmetricCscMatrix(true), true, SuiteSparseOrdering.AMD))
            //{
            //    watch.Stop();
            //    Console.WriteLine($"{unorderedName} ordering -> factorization (with hidden AMD): Non zeros = {factor.NumNonZerosUpper}"
            //        + $" , time for factorization = {watch.ElapsedMilliseconds} ms");

            //}

            //if (printEnumerations)
            //{
            //    Console.WriteLine("\nPermutation: ");
            //    PrintList(permutation);
            //    Console.WriteLine("\nBefore reordering: ");
            //    unorderedDofs.WriteToConsole();
            //    Console.WriteLine("\nAfter reordering: ");
            //    reorderedDofs.WriteToConsole();
            //}
        }

        private static Model2D CreateDcb3x1()
        {
            var benchmark = new DCB3x1(20);
            Model2D model = benchmark.CreateModel();
            benchmark.HandleEnrichment(model);
            return model;
        }

        private static Model2D CreateDcbGRACM()
        {
            double growthLength = 0.3;
            double fineElementSize = 0.045;
            IMeshProvider meshProvider = new DCBRefinedMeshProvider(fineElementSize, 10 * fineElementSize);
            //IMeshProvider meshProvider = new GmshMeshProvider(@"C: \Users\Serafeim\Desktop\GMSH\dcb.msh");
            var builder = new DCB.Builder(growthLength, meshProvider);
            builder.UseLSM = true;
            //TODO: fix a bug that happens when the crack has almost reached the boundary, is inside but no tip elements are 
            //      found. It happens at iteration 10.
            builder.MaxIterations = 10;
            builder.KnownPropagation = null; // TODO: enter the fixed propagator here, perhaps by solving the benchmark once.


            var benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();

            benchmark.Crack.UpdateGeometry(0.0475718018434623, 0.3);
            benchmark.Crack.UpdateEnrichments();
            benchmark.Crack.UpdateGeometry(-0.0694073863911697, 0.3);
            benchmark.Crack.UpdateEnrichments();
            benchmark.Crack.UpdateGeometry(-0.0963237968537083, 0.3);
            benchmark.Crack.UpdateEnrichments();
            benchmark.Crack.UpdateGeometry(-0.118290637917496, 0.3);
            benchmark.Crack.UpdateEnrichments();

            return benchmark.Model;
        }

        private static void PrintList(IReadOnlyList<int> list)
        {
            for (int i = 0; i < list.Count; ++i) Console.Write(list[i] + " ");
            Console.WriteLine();
        }
    }
}
