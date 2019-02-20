using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class SubdomainConstrainedDofOrderingGeneral : ISubdomainConstrainedDofOrdering
    {
        public SubdomainConstrainedDofOrderingGeneral(int numConstrainedDofs, DofTable subdomainConstrainedDofs)
        {
            this.NumConstrainedDofs = numConstrainedDofs;
            this.ConstrainedDofs = subdomainConstrainedDofs;
        }

        public DofTable ConstrainedDofs { get; }
        public int NumConstrainedDofs { get; }
        
        public (int[] elementDofIndices, int[] subdomainDofIndices) MapConstrainedDofsElementToSubdomain(IElement_v2 element)
        {
            IList<INode> elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            IList<IList<DOFType>> elementDofs = element.ElementType.DofEnumerator.GetDOFTypes(element);

            // Count the dof superset (free and constrained) to allocate enough memory and avoid resizing
            int allElementDofs = 0;
            for (int i = 0; i < elementNodes.Count; ++i) allElementDofs += elementDofs[i].Count;
            var elementDofIndices = new List<int>(allElementDofs);
            var subdomainDofIndices = new List<int>(allElementDofs);

            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isConstrained = ConstrainedDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx],
                        out int subdomainDofIdx);
                    if (isConstrained)
                    {
                        elementDofIndices.Add(elementDofIdx);
                        subdomainDofIndices.Add(subdomainDofIdx);
                    }
                    ++elementDofIdx; // This must be incremented for free dofs as well
                }
            }
            return (elementDofIndices.ToArray(), subdomainDofIndices.ToArray());
        }
    }
}
