using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Tests.Khoei;

namespace ISAAR.MSolve.XFEM.Tests.DofOrdering
{
    static class TestEnumerator
    {
        public static void Run()
        {
            var dcb = new DCB3x1(20);
            Model2D model = dcb.CreateModel();
            dcb.HandleEnrichment(model);
            TestReorderDeepCopy(model);
            PrintElementDofs(model);
        }

        /// <summary>
        /// The reordering just increases each number by 1. This method should only change the copied DofOrderer
        /// </summary>
        /// <param name="model"></param>
        /// <param name="solver"></param>
        private static void TestReorderDeepCopy(Model2D model)
        {
            var enumerator = InterleavedDofOrderer.Create(model);
            var permutationOldToNew = new SortedDictionary<int, int>();
            foreach (XNode2D node in model.Nodes)
            {
                // Standard X dof
                try
                {
                    int idxOldX = enumerator.GetStandardDofOf(node, DisplacementDof.X);
                    permutationOldToNew.Add(idxOldX, idxOldX + 1);
                }
                catch (KeyNotFoundException)
                { }


                // Standard Y dof
                try
                {
                    int idxOldY = enumerator.GetStandardDofOf(node, DisplacementDof.Y);
                    permutationOldToNew.Add(idxOldY, idxOldY + 1);
                }
                catch (KeyNotFoundException)
                { }

                foreach (var enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (var dof in enrichment.Dofs)
                    {
                        int idxOld = enumerator.GetEnrichedDofOf(node, dof);
                        permutationOldToNew.Add(idxOld, idxOld + 1);
                    }
                }
            }
            
            Console.WriteLine("\n------------------- Renumbering ---------------------");
            Console.WriteLine("Before renumbering:");
            enumerator.WriteToConsole();
            
            IDofOrderer enumeratorCopy = enumerator.DeepCopy();
            enumeratorCopy.ReorderUnconstrainedDofs(permutationOldToNew.Values.ToArray(), true);
            Console.WriteLine("\nAfter renumbering original enumerator:");
            enumerator.WriteToConsole(); // They should not have changed
            Console.WriteLine("\nAfter renumbering copied and reordered enumerator:");
            enumeratorCopy.WriteToConsole();
        }

        public static void PrintElementDofs(Model2D model)
        {
            var enumerator = SeparateDofOrderer.Create(model);
            int counter = 0;
            foreach (var element in model.Elements)
            {
                Console.WriteLine($"Element {counter++}:");
                Console.Write("Free standard dofs: ");
                PrintList(enumerator.GetStandardDofsOf(element));
                Console.Write("Constrained standard dofs: ");
                PrintList(enumerator.GetConstrainedDofsOf(element));
                Console.Write("Free enriched dofs: ");
                PrintList(enumerator.GetEnrichedDofsOf(element));
                Console.WriteLine();
            }
        }

        private static void PrintList(IReadOnlyList<int> list)
        {
            for (int i = 0; i < list.Count; ++i) Console.Write(list[i] + " ");
            Console.WriteLine();
        }
    }
}
