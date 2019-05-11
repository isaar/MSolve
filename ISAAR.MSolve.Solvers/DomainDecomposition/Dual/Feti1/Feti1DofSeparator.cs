using System.Collections.Generic;
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

        internal void SeparateDofs(IStructuralModel model)
        {
            GatherGlobalBoundaryDofs(model.Nodes, model.GlobalDofOrdering);
            SeparateBoundaryInternalDofs(model);
        }

        private void SeparateBoundaryInternalDofs(IStructuralModel model)
        {
            InternalDofIndices = new Dictionary<int, int[]>();
            BoundaryDofIndices = new Dictionary<int, int[]>();
            BoundaryDofs = new Dictionary<int, (INode node, IDofType dofType)[]>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                (int[] internalDofIndices, int[] boundaryDofIndices, (INode node, IDofType dofType)[] boundaryDofConnectivities)
                    = base.SeparateBoundaryInternalDofs(subdomain.Nodes, subdomain.FreeDofOrdering.FreeDofs);

                InternalDofIndices.Add(subdomain.ID, internalDofIndices);
                BoundaryDofIndices.Add(subdomain.ID, boundaryDofIndices);
                BoundaryDofs.Add(subdomain.ID, boundaryDofConnectivities);
            }
        }
    }
}
