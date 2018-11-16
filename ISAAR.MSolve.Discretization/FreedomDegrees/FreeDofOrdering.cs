using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class FreeDofOrdering: IDofOrdering
    {
        public FreeDofOrdering(int numFreeDofs, DofTable freeDofs)
        {
            this.NumFreeDofs = numFreeDofs;
            this.FreeDofs = freeDofs;
        }

        public DofTable FreeDofs { get; }
        public int NumFreeDofs { get; }

        public IVector ExtractVectorElementFromGlobal(IElement element, IVectorView globalVector)
        {
            IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
            IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypes(element);
            int numElementDofs = 0;
            for (int nodeIdx = 0; nodeIdx < elementDofs.Count; ++nodeIdx) numElementDofs += elementDofs[nodeIdx].Count;

            double[] elementVector = new double[numElementDofs];
            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx],
                        out int globalDofIdx);
                    if (isFree) elementVector[elementDofIdx] = globalVector[globalDofIdx];
                    // Else, the quantity of interest is 0.0 at all constrained dofs.
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
            return Vector.CreateFromArray(elementVector);
        }

        public IReadOnlyDictionary<int, int> MapFreeDofsElementToGlobal(IElement element)
        {
            IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
            IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypes(element);

            var dofMap = new Dictionary<int, int>();
            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx],
                        out int globalDofIdx);
                    if (isFree) dofMap[elementDofIdx] = globalDofIdx;
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
            return dofMap;
        }
    }
}
