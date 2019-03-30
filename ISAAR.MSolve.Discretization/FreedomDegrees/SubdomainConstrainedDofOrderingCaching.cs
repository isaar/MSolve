using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: This should be thread safe.
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class SubdomainConstrainedDofOrderingCaching : ISubdomainConstrainedDofOrdering
    {
        private readonly Dictionary<IElement_v2, (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices)> 
            elementDofsCache = new Dictionary<IElement_v2, (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices)>();

        public SubdomainConstrainedDofOrderingCaching(int numConstrainedDofs, DofTable subdomainConstrainedDofs)
        {
            this.NumConstrainedDofs = numConstrainedDofs;
            this.ConstrainedDofs = subdomainConstrainedDofs;
        }

        public DofTable ConstrainedDofs { get; }
        public int NumConstrainedDofs { get; }

        public (int[] elementDofIndices, int[] subdomainDofIndices) MapConstrainedDofsElementToSubdomain(IElement_v2 element)
        {
            (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices) = GetElementData(element);
            return (elementDofIndices, subdomainDofIndices);
        }

        private (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices) GetElementData(IElement_v2 element)
        {
            bool isStored = elementDofsCache.TryGetValue(element, out (int, int[], int[]) elementData);
            if (isStored) return elementData;
            else
            {
                elementData = ProcessElement(element);
                elementDofsCache.Add(element, elementData);
                return elementData;
            }
        }

        private (int numAllElementDofs, int[] elementDofIndices, int[] subdomainDofIndices) ProcessElement(IElement_v2 element)
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
            return (allElementDofs, elementDofIndices.ToArray(), subdomainDofIndices.ToArray());
        }
    }
}
