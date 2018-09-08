using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class EnrichedOrderingAmd : IEnrichedOrdering
    {
        public void ReorderEnrichedDofs(XSubdomain2D subdomain)
        {
            var orderingAlgorithm = new OrderingAmd();

            int order = subdomain.DofOrderer.NumEnrichedDofs;
            var pattern = SparsityPatternSymmetric.CreateEmpty(order);
            // Could build the sparsity pattern during Dof enumeration?
            foreach (var element in subdomain.Elements)
            {
                var enrichedDofs = subdomain.DofOrderer.GetSubdomainEnrichedDofsOf(element);
                pattern.ConnectIndices(enrichedDofs, false);
            }
            (int[] permutation, ReorderingStatistics stats) = orderingAlgorithm.FindPermutation(pattern);

            subdomain.DofOrderer.ReorderSubdomainDofs(permutation, false);
        }
    }
}
