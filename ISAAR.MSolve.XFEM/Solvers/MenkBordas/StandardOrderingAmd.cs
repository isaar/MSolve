using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class StandardOrderingAmd : IStandardOrdering
    {
        private readonly Model2D model;

        public StandardOrderingAmd(Model2D model)
        {
            this.model = model;
        }

        public void ReorderStandardDofs(XClusterDofOrderer stdDofOrderer)
        {
            var orderingAlgorithm = new OrderingAmd();

            int order = stdDofOrderer.NumStandardDofs;
            var pattern = SparsityPatternSymmetric.CreateEmpty(order);
            // Could build the sparsity pattern during Dof enumeration?
            foreach (var element in model.Elements)
            {
                var standardDofs = stdDofOrderer.GetStandardDofsOf(element);
                pattern.ConnectIndices(standardDofs, false);
            }
            (int[] permutation, ReorderingStatistics stats) = orderingAlgorithm.FindPermutation(pattern);

            stdDofOrderer.ReorderStandardDofs(permutation, false);
        }
    }
}
