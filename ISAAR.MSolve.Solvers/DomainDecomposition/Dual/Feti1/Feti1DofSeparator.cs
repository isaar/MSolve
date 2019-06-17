using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices;

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

        public void DefineGlobalBoundaryDofs(IStructuralModel model)
        {
            base.DefineGlobalBoundaryDofs(model.Nodes, model.GlobalDofOrdering);
        }

        public void SeparateBoundaryInternalDofs(ISubdomain subdomain)
        {
            int s = subdomain.ID;
            (int[] internalDofIndices, int[] boundaryDofIndices, (INode node, IDofType dofType)[] boundaryDofConnectivities)
                = DofSeparatorBase.SeparateBoundaryInternalDofs(subdomain.Nodes, subdomain.FreeDofOrdering.FreeDofs);

            InternalDofIndices[s] = internalDofIndices;
            BoundaryDofIndices[s] = boundaryDofIndices;
            BoundaryDofs[s] = boundaryDofConnectivities;
        }
    }
}
