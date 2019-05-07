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
        //TODO: BoundaryDofConnectivities is useful only for heterogeneous problems and its creation might be costly.
        //      Needs benchmarking and possibly a way to avoid it in homogeneous problems
        internal Dictionary<int, (INode node, IDofType dofType)[]> BoundaryDofConnectivities { get; private set; }
        internal Dictionary<int, int[]> BoundaryDofIndices { get; private set; }
        internal Dictionary<int, int[]> BoundaryDofMultiplicities { get; private set; }
        internal Dictionary<int, int[]> InternalDofIndices { get; private set; }

        internal void SeparateDofs(IStructuralModel model)
        {
            base.GatherDualDofs(model.Nodes, model.GlobalDofOrdering);
            SeparateBoundaryInternalDofs(model);
        }

        private void SeparateBoundaryInternalDofs(IStructuralModel model)
        {
            InternalDofIndices = new Dictionary<int, int[]>();
            BoundaryDofIndices = new Dictionary<int, int[]>();
            BoundaryDofMultiplicities = new Dictionary<int, int[]>();
            BoundaryDofConnectivities = new Dictionary<int, (INode node, IDofType dofType)[]>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                (int[] internalDofIndices, int[] boundaryDofIndices, int[] boundaryDofMultiplicities, 
                    (INode node, IDofType dofType)[] boundaryDofConnectivities) = 
                    base.SeparateBoundaryInternalDofs(subdomain.Nodes, subdomain.FreeDofOrdering.FreeDofs);

                InternalDofIndices.Add(subdomain.ID, internalDofIndices);
                BoundaryDofIndices.Add(subdomain.ID, boundaryDofIndices);
                BoundaryDofMultiplicities.Add(subdomain.ID, boundaryDofMultiplicities);
                BoundaryDofConnectivities.Add(subdomain.ID, boundaryDofConnectivities);
            }
        }
    }
}
