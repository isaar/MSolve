using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: This is dramatically slower than NodeMajorDofOrderingStrategy. It must be made faster. Also this + NodeMajorReordering() 
//      must be at least as fast as NodeMajorDofOrderingStrategy. Then the solvers should have simple + reordering as defaults.
namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// WARNING: DO NOT USE THIS YET. IT IS VERY SLOW.
    /// Free dofs are assigned global / subdomain indices according to the order they are first encountered. 
    /// Constrained dofs are ignored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SimpleDofOrderingStrategy //: IFreeDofOrderingStrategy
    {
        public (int numGlobalFreeDofs, DofTable globalFreeDofs) OrderGlobalDofs(IStructuralModel model)
            => OrderFreeDofsOfElementSet(model.Elements, model.Constraints);

        public (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) OrderSubdomainDofs(ISubdomain subdomain)
            => OrderFreeDofsOfElementSet(subdomain.Elements, subdomain.Constraints);

        private static (int numFreeDofs, DofTable freeDofs) OrderFreeDofsOfElementSet(IEnumerable<IElement> elements,
            Table<INode, IDofType, double> constraints)
        {
            var freeDofs = new DofTable();
            int dofCounter = 0;
            foreach (IElement element in elements)
            {
                //IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element); //this is wrong
                IReadOnlyList<INode> elementNodes = element.Nodes;
                IReadOnlyList<IReadOnlyList<IDofType>> elementDofs = element.ElementType.DofEnumerator.GetDofTypesForDofEnumeration(element);
                for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
                {
                    bool isNodeConstrained = constraints.TryGetDataOfRow(elementNodes[nodeIdx],
                        out IReadOnlyDictionary<IDofType, double> constraintsOfNode);
                    for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                    {
                        IDofType dofType = elementDofs[nodeIdx][dofIdx];
                        bool isDofConstrained = isNodeConstrained ? constraintsOfNode.ContainsKey(dofType) : false;
                        if (!isDofConstrained)
                        {
                            bool isNewDof = freeDofs.TryAdd(elementNodes[nodeIdx], dofType, dofCounter);
                            if (isNewDof) ++dofCounter;
                        }
                    }
                }
            }
            return (dofCounter, freeDofs);
        }
    }
}
