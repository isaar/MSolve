using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.XFEM.Entities;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: Find out why it reports that dense rows are moved to the end for fillet benchmark
//TODO: compare with unconstrained AMD. From some small tests CAMD is generally better, which doesn't make sense.
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class EnrichedOrderingCamd: IEnrichedOrdering
    {
        public void ReorderEnrichedDofs(XSubdomain2D subdomain)
        {

            int order = subdomain.DofOrderer.NumEnrichedDofs;
            var pattern = SparsityPatternSymmetric.CreateEmpty(order);
            
            // Sparsity pattern of U
            foreach (var element in subdomain.Elements)
            {
                var enrichedDofs = subdomain.DofOrderer.GetSubdomainEnrichedDofsOf(element);
                pattern.ConnectIndices(enrichedDofs, false);
            }

            
            int[] boundaryDofs = FindEnrichedBoundaryDofs(subdomain);
            if (boundaryDofs.Length == 0)
            {
                var amd = new OrderingAmd();
                (int[] permutation, ReorderingStatistics stats) = amd.FindPermutation(pattern);
                subdomain.DofOrderer.ReorderSubdomainDofs(permutation, false);
            }
            else // Constraints: move boundary dofs to the end
            {
                var camd = new OrderingCamd(-1, true); // Make sure dense rows are not moved
                var constraints = new int[order]; // internal dofs have ordinal = 0, to end up first
                foreach (var dof in boundaryDofs) constraints[dof] = 1; // boundary dofs have ordinal = 1

                (int[] permutation, ReorderingStatistics stats) = camd.FindPermutation(pattern, constraints);
                if (stats.NumMovedDenseRows != 0) //TODO: create a custom exception class
                {
                    //throw new Exception("In order for enriched boundary dofs to be moved to the end, no dense rows must not be"
                    //    + $" moved there themselves. However {stats.NumMovedDenseRows} dense rows were moved to the end"); 


                    //TODO: Find out why it reports that dense rows are moved to the end for fillet benchmark
                }
                subdomain.DofOrderer.ReorderSubdomainDofs(permutation, false);

                #region debug
                //CheckBoundaryDofsAreAtTheEnd(order, boundaryDofs, permutation);
                //Console.Write("permutation (new-to-old) = ");
                //for (int i = 0; i < order; ++i) Console.Write(permutation[i] + " ");
                //Console.WriteLine();
                #endregion
            }
        }

        private void CheckBoundaryDofsAreAtTheEnd(int order, int[] boundaryDofs, int[] permutation)
        {
            var lastDofs = new HashSet<int>(); // Contains the original indices of these dofs
            int lastDofsStart = order - boundaryDofs.Length;
            for (int i = lastDofsStart; i < order; ++i) lastDofs.Add(permutation[i]);
            foreach (var dof in boundaryDofs) //TODO: create a custom exception class
            {
                if (!lastDofs.Contains(dof)) throw new Exception($"Boundary dof {dof} was not placed at the end");
            }
        }

        private static int[] FindEnrichedBoundaryDofs(XSubdomain2D subdomain)
        {
            int numBoundaryDofs = 0;
            foreach (var nodeDofsPair in subdomain.DofOrderer.BoundaryDofs) numBoundaryDofs += nodeDofsPair.Value.Count;

            var boundaryDofs = new int[numBoundaryDofs];
            int idx = 0;
            foreach (var nodeDofsPair in subdomain.DofOrderer.BoundaryDofs)
            {
                XNode2D node = nodeDofsPair.Key;
                foreach (var dof in nodeDofsPair.Value)
                {
                    boundaryDofs[idx++] = subdomain.DofOrderer.GetSubdomainEnrichedDofOf(node, dof);
                }
            }

            #region debug
            //Console.Write($"Subdomain {subdomain.ID}: boundary enriched dofs: ");
            //if (numBoundaryDofs == 0) Console.Write("none");
            //for (int i = 0; i < numBoundaryDofs; ++i) Console.Write(boundaryDofs[i] + " ");
            //Console.WriteLine();
            #endregion

            return boundaryDofs;
        }
    }
}
