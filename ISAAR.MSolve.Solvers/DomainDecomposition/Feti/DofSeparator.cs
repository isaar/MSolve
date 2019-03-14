using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Needs a proper name. This probably cannot be incorporated in the ISubdomainDofOrdering, as the intent is different and
//      depending on the DD method the dof categories may be different (e.g. FETI-1: internal/boundary, 
//      FETI-DP: corner/boundary/remainder)
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public class DofSeparator
    {
        //TODO: each subdomain should have an object that manages these arrays.
        private Dictionary<int, int[]> boundaryDofs;
        private Dictionary<int, int[]> boundaryDofsMultiplicity;
        private Dictionary<int, int[]> internalDofs;

        internal void SeparateBoundaryInternalDofs(Dictionary<int, ISubdomain_v2> subdomains)
        {
            boundaryDofs = new Dictionary<int, int[]>();
            boundaryDofsMultiplicity = new Dictionary<int, int[]>();
            internalDofs = new Dictionary<int, int[]>();
            foreach (ISubdomain_v2 subdomain in subdomains.Values)
            {
                var boundaryDofsOfSubdomain = new SortedDictionary<int, int>(); // key = dofIdx, value = multiplicity
                var internalDofsOfSubdomain = new SortedSet<int>();
                foreach (INode node in subdomain.Nodes)
                {
                    IEnumerable<int> nodalDofs = subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node);
                    int nodeMultiplicity = node.SubdomainsDictionary.Count;
                    if (nodeMultiplicity > 1) // boundary node
                    {
                        foreach (int dof in nodalDofs) boundaryDofsOfSubdomain.Add(dof, nodeMultiplicity);
                    }
                    else
                    {
                        foreach (int dof in nodalDofs) internalDofsOfSubdomain.Add(dof);
                    }
                }
                boundaryDofs.Add(subdomain.ID, boundaryDofsOfSubdomain.Keys.ToArray());
                boundaryDofsMultiplicity.Add(subdomain.ID, boundaryDofsOfSubdomain.Values.ToArray()); // sorted the same as Keys
                internalDofs.Add(subdomain.ID, internalDofsOfSubdomain.ToArray());
            }
        }

        internal Dictionary<int, int[]> BoundaryDofs => boundaryDofs;
        internal Dictionary<int, int[]> BoundaryDofsMultiplicity => boundaryDofsMultiplicity;
        internal Dictionary<int, int[]> InternalDofs => internalDofs;
    }
}
