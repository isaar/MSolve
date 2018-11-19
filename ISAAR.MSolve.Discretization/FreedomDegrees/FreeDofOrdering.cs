using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class FreeDofOrdering: IDofOrdering
    {
        public FreeDofOrdering(int numFreeDofs, DofTable freeDofs, Dictionary<int, Dictionary<DOFType, int>> globalFreeDofs)
        {
            this.NumFreeDofs = numFreeDofs;
            this.FreeDofs = freeDofs;

            FreeDofMapSubdomainToGlobal = new int[numFreeDofs];
            foreach ((INode node, DOFType dofType, int subdomainDofIdx) in freeDofs)
            {
                FreeDofMapSubdomainToGlobal[subdomainDofIdx] = globalFreeDofs[node.ID][dofType];
            }
        }

        public int[] FreeDofMapSubdomainToGlobal { get; }
        public DofTable FreeDofs { get; }
        public int NumFreeDofs { get; }

        public void AddVectorElementToSubdomain(IElement element, IVectorView elementVector, IVector subdomainVector)
        {
            IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
            IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypes(element);

            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx],
                        out int subdomainDofIdx);
                    if (isFree)
                    {
                        subdomainVector.Set(subdomainDofIdx, subdomainVector[subdomainDofIdx] + elementVector[elementDofIdx]);
                    }
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
        }

        public void AddVectorSubdomainToGlobal(ISubdomain_v2 subdomain, IVectorView subdomainVector, IVector globalVector)
        {
            for (int i = 0; i < NumFreeDofs; ++i)
            {
                int globalIdx = FreeDofMapSubdomainToGlobal[i];
                globalVector.Set(globalIdx, subdomainVector[i] + globalVector[globalIdx]);
                //TODO: add a Vector.SetSubvector and Vector.AddSubvector for incontiguous entries
            }

            //foreach (INode node in subdomain.Nodes)
            //{
            //    //Dictionary<DOFType, int> dofTypes = NodalDOFsDictionary[nodeID];
            //    foreach (DOFType dofType in FreeDofs.GetColumnsOfRow(node))
            //    {
            //        bool isFreeDof = FreeDofs.TryGetValue(node, dofType, out int subdomainDofIdx);
            //        if (isFreeDof)
            //        {
            //            //int globalDofIdx = subdomain.GlobalNodalDOFsDictionary[node.ID][dofType];
            //            //Debug.Assert(globalDofIdx > -1);

            //            int globalDofIdx = subdomain.GlobalFreeDofs[node, dofType];

            //            //TODO: add a Vector.SetSubvector and Vector.AddSubvector for incontiguous entries
            //            globalVector.Set(globalDofIdx, globalVector[globalDofIdx] + subdomainVector[subdomainDofIdx]);
            //        }
            //    }
            //}
        }

        public double[] ExtractVectorElementFromSubdomain(IElement element, IVectorView subdomainVector)
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
                        out int subdomainDofIdx);
                    if (isFree) elementVector[elementDofIdx] = subdomainVector[subdomainDofIdx];
                    // Else, the quantity of interest is 0.0 at all constrained dofs.
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
            return elementVector;
        }

        //TODO: this needs refactoring: DofTable for the global dofs, perhaps remove the global dofs from the subdomain, ...
        public void ExtractVectorSubdomainFromGlobal(ISubdomain_v2 subdomain, IVectorView globalVector, IVector subdomainVector)
        {
            for (int i = 0; i < NumFreeDofs; ++i)
            {
                //TODO: add a Vector.SetSubvector and Vector.AddSubvector for incontiguous entries
                subdomainVector.Set(i, globalVector[FreeDofMapSubdomainToGlobal[i]]);
            }

            //foreach (INode node in subdomain.Nodes)
            //{
            //    foreach (DOFType dofType in FreeDofs.GetColumnsOfRow(node))
            //    {
            //        bool isFreeDof = FreeDofs.TryGetValue(node, dofType, out int subdomainDodIdx);

            //        if (isFreeDof)
            //        {
            //            //int globalDofIdx = subdomain.GlobalNodalDOFsDictionary[node.ID][dofType];
            //            //Debug.Assert(globalDofIdx > -1);

            //            int globalDofIdx = subdomain.GlobalFreeDofs[node, dofType];


            //            //TODO: add a Vector.SetSubvector and Vector.AddSubvector for incontiguous entries
            //            subdomainVector.Set(subdomainDodIdx, globalVector[globalDofIdx]);
            //        }
            //    }
            //}
        }

        public IReadOnlyDictionary<int, int> MapFreeDofsElementToSubdomain(IElement element)
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
                        out int subdomainDofIdx);
                    if (isFree) dofMap[elementDofIdx] = subdomainDofIdx;
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
            return dofMap;
        }
    }
}
