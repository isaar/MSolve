using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Needs a proper name. This probably cannot be incorporated in the ISubdomainDofOrdering, as the intent is different and
//      depending on the DD method the dof categories may be different (e.g. FETI-1: internal/boundary, 
//      FETI-DP: corner/boundary/remainder)
//TODO: Not sure about having the indexing data of all subdomains in the same class.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1
{
    public class Feti1DofSeparator : DofSeparatorBase
    {
        public override Dictionary<int, int[]> BoundaryDofIndices { get; protected set; }
        public override Dictionary<int, (INode node, IDofType dofType)[]> BoundaryDofs { get; protected set; }
        public override Dictionary<int, int[]> InternalDofIndices { get; protected set; }

        public Feti1DofSeparator()
        {
            InternalDofIndices = new Dictionary<int, int[]>();
            BoundaryDofIndices = new Dictionary<int, int[]>();
            BoundaryDofs = new Dictionary<int, (INode node, IDofType dofType)[]>();
        }

        internal void SeparateDofs(IStructuralModel model)
        {
            GatherGlobalBoundaryDofs(model.Nodes, model.GlobalDofOrdering);
            SeparateBoundaryInternalDofs(model);
        }

        private void SeparateBoundaryInternalDofs(IStructuralModel model)
        {
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                if (!subdomain.ConnectivityModified) continue;

                int s = subdomain.ID;
                Debug.WriteLine($"Separating boundary-internal dofs of subdomain {s}");
                (int[] internalDofIndices, int[] boundaryDofIndices, (INode node, IDofType dofType)[] boundaryDofConnectivities)
                    = base.SeparateBoundaryInternalDofs(subdomain.Nodes, subdomain.FreeDofOrdering.FreeDofs);

                InternalDofIndices[s] = internalDofIndices;
                BoundaryDofIndices[s] = boundaryDofIndices;
                BoundaryDofs[s] = boundaryDofConnectivities;
            }
        }
    }
}
