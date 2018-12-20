using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;

//TODO: The static factory methods are fine for now (and as defaults), but if the algorithm classes need options, then the user 
//      must instantiate them. There should probably be an IAmdAlgorithm interface. Or this class should not be named AMD.
namespace ISAAR.MSolve.Solvers.Ordering.Reordering
{
    /// <summary>
    /// Reorders the unconstrained freedom degrees according to the fill-reducing permutation calculated by the Approximate 
    /// Minimum Degree algorithm. Note that the pattern of the sparse matrix, i.e. the positions of its non-zero entries, must 
    /// be constructed and then passed to AMD. These might be costly operations and AMD might fail.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class AmdReordering : IDofReorderingStrategy
    {
        private readonly IReorderingAlgorithm amd;

        private AmdReordering(IReorderingAlgorithm amd)
        {
            this.amd = amd;
        }

        public static AmdReordering CreateWithSuiteSparseAmd() => new AmdReordering(new OrderingAmdSuiteSparse());
        public static AmdReordering CreateWithCSparseAmd() => new AmdReordering(new OrderingAmdCSparseNet());

        public void ReorderDofs(ISubdomain_v2 subdomain, ISubdomainFreeDofOrdering originalOrdering)
        {
            originalOrdering.Reorder(amd, subdomain);

            //TODO: Delete the next. That code has been copied to the ISubdomainFreeDofOrdering implementations. This class
            //      works for any reordering algorithm, not only AMD.
            //var pattern = SparsityPatternSymmetric.CreateEmpty(originalOrdering.NumFreeDofs);
            //foreach (var element in subdomain.Elements)
            //{
            //    (int[] elementDofIndices, int[] subdomainDofIndices) = originalOrdering.MapFreeDofsElementToSubdomain(element);

            //    //TODO: ISubdomainFreeDofOrdering could perhaps return whether the subdomainDofIndices are sorted or not.
            //    pattern.ConnectIndices(subdomainDofIndices, false);
            //}
            //(int[] permutation, bool oldToNew) = amd.FindPermutation(pattern);
        }
    }
}
