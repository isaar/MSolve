using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: This should be thread safe.
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class SubdomainFreeDofOrderingCaching : ISubdomainFreeDofOrdering
    {
        private readonly Dictionary<IElement, (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices)> 
            elementDofsCache = new Dictionary<IElement, (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices)>();

        public SubdomainFreeDofOrderingCaching(int numFreeDofs, DofTable subdomainFreeDofs)
        {
            this.NumFreeDofs = numFreeDofs;
            this.FreeDofs = subdomainFreeDofs;
        }

        public DofTable FreeDofs { get; }
        public int NumFreeDofs { get; }

        public void AddVectorElementToSubdomain(IElement element, IVectorView elementVector, IVector subdomainVector)
        {
            (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices) = GetElementData(element);
            subdomainVector.AddNonContiguouslyFrom(subdomainDofIndices, elementVector, elementDofIndices);
        }

        public int CountElementDofs(IElement element)
        {
            (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices) = GetElementData(element);
            return numAllDofs;
        }

        public void ExtractVectorElementFromSubdomain(IElement element, IVectorView subdomainVector, IVector elementVector)
        {
            (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices) = GetElementData(element);

            elementVector.CopyNonContiguouslyFrom(elementDofIndices, subdomainVector, subdomainDofIndices);

            //double[] elementVector = new double[numAllDofs];
            //for (int i = 0; i < elementDofIndices.Length; ++i)
            //{
            //    int elementDofIdx = elementDofIndices[i];
            //    int subdomainDofIdx = subdomainDofIndices[i];
            //    elementVector[elementDofIdx] = subdomainVector[subdomainDofIdx];
            //}
            //return elementVector;
        }

        public (int[] elementDofIndices, int[] subdomainDofIndices) MapFreeDofsElementToSubdomain(IElement element)
        {
            (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices) = GetElementData(element);
            return (elementDofIndices, subdomainDofIndices);
        }

        public void Reorder(IReorderingAlgorithm reorderingAlgorithm, ISubdomain_v2 subdomain)
        {
            elementDofsCache.Clear();
            var pattern = SparsityPatternSymmetric.CreateEmpty(NumFreeDofs);
            foreach (var element in subdomain.Elements)
            {
                // Do not cache anything at this point
                (int numAllElementDofs, int[] elementDofIndices, int[] subdomainDofIndices) = ProcessElement(element);

                //TODO: ISubdomainFreeDofOrdering could perhaps return whether the subdomainDofIndices are sorted or not.
                pattern.ConnectIndices(subdomainDofIndices, false);
            }
            (int[] permutation, bool oldToNew) = reorderingAlgorithm.FindPermutation(pattern);
            FreeDofs.Reorder(permutation, oldToNew);
        }

        public void ReorderNodeMajor(IReadOnlyList<INode> sortedNodes)
        {
            elementDofsCache.Clear();
            FreeDofs.ReorderNodeMajor(sortedNodes);
        }

        private (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices) GetElementData(IElement element)
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

        private (int numAllElementDofs, int[] elementDofIndices, int[] subdomainDofIndices) ProcessElement(IElement element)
        {
            IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
            IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypes(element);

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
                    bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx],
                        out int subdomainDofIdx);
                    if (isFree)
                    {
                        elementDofIndices.Add(elementDofIdx);
                        subdomainDofIndices.Add(subdomainDofIdx);
                    }
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
            return (allElementDofs, elementDofIndices.ToArray(), subdomainDofIndices.ToArray());
        }
    }
}
