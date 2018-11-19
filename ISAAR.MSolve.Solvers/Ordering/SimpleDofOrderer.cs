using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Free dofs are assigned global (actually subdomain) indices according to the order they are first encountered. 
    /// Constrained dofs are ignored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SimpleDofOrderer: IDofOrderer
    {
        public IDofOrdering OrderDofs(ISubdomain_v2 subdomain)
        {
            var freeDofs = new DofTable();
            int dofCounter = 0;
            foreach (IElement element in subdomain.ΙElementsDictionary.Values)
            {
                //IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element); //this is wrong
                IList<INode> elementNodes = element.INodes;
                IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element);
                for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
                {
                    bool isNodeConstrained = subdomain.Constraints.TryGetDataOfRow(elementNodes[nodeIdx],
                        out IReadOnlyDictionary<DOFType, double> constraintsOfNode);
                    for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                    {
                        DOFType dofType = elementDofs[nodeIdx][dofIdx];
                        bool isDofConstrained = isNodeConstrained ? constraintsOfNode.ContainsKey(dofType) : false;
                        if (!isDofConstrained)
                        {
                            bool isNewDof = freeDofs.TryAdd(elementNodes[nodeIdx], dofType, dofCounter);
                            if (isNewDof) ++dofCounter;
                        }
                    }
                }
            }
            return new FreeDofOrdering(dofCounter, freeDofs);
        }
    }
}
