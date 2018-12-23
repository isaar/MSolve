using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering.Reordering
{
    /// <summary>
    /// Reorders the unconstrained freedom degrees of a subdomain, such that they are sorted in a node major fashion: The dofs 
    /// of the first node are numbered, then the dofs of the second node, etc.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class NodeMajorReordering : IDofReorderingStrategy
    {
        public void ReorderDofs(ISubdomain_v2 subdomain, ISubdomainFreeDofOrdering originalOrdering)
            => originalOrdering.ReorderNodeMajor(subdomain.Nodes);
    }
}
