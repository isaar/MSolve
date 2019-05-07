using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Should this be abstract or not suffixed with "Base".
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual
{
    public class DofSeparatorBase : IDofSeparator
    {
        public Dictionary<INode, IDofType[]> DualDofs { get; private set; }

        protected void GatherDualDofs(IEnumerable<INode> allNodes, IGlobalFreeDofOrdering globalDofOrdering)
        {
            DualDofs = new Dictionary<INode, IDofType[]>();

            //TODO: model.Nodes probably doesn't work if there are embedded nodes. It is time to isolate the embedded nodes. Or I could use the GlobalDofOrdering.
            foreach (INode node in allNodes)
            {
                int nodeMultiplicity = node.SubdomainsDictionary.Count;
                if (nodeMultiplicity > 1)
                {
                    // Access the free dofs only. Does this also filter out embedded dofs?
                    IDofType[] dofsOfNode = globalDofOrdering.GlobalFreeDofs.GetColumnsOfRow(node).ToArray(); //TODO: interacting with the table can be optimized.

                    // If all dofs of this node are constrained, then it is not considered boundary.
                    if (dofsOfNode.Length == 0) continue;
                    else DualDofs[node] = dofsOfNode;
                }
            }
        }

        protected (int[] internalDofIndices, int[] boundaryDofIndices, int[] boundaryDofMultiplicities, 
            (INode node, IDofType dofType)[] boundaryDofConnectivities) 
            SeparateBoundaryInternalDofs(IEnumerable<INode> nodes, DofTable freeDofs)
        {
            var boundaryMultiplicities = new SortedDictionary<int, int>(); // key = dofIdx, value = multiplicity
            var boundaryConnectivities = new SortedDictionary<int, (INode node, IDofType dofType)>();
            var internalDofs = new SortedSet<int>(); //TODO: Why set instead of List?
            foreach (INode node in nodes)
            {
                int nodeMultiplicity = node.SubdomainsDictionary.Count;
                if (nodeMultiplicity > 1) // boundary node
                {
                    bool isNodeFree = freeDofs.TryGetDataOfRow(node, 
                        out IReadOnlyDictionary<IDofType, int> dofTypesIndices); // This avoids embedded and constrained dofs
                    if (isNodeFree)
                    {
                        foreach (var dofTypeIdxPair in dofTypesIndices)
                        {
                            boundaryMultiplicities.Add(dofTypeIdxPair.Value, nodeMultiplicity);
                            boundaryConnectivities.Add(dofTypeIdxPair.Value, (node, dofTypeIdxPair.Key));
                        }
                    }
                }
                else // internal nodes
                {
                    foreach (int dof in freeDofs.GetValuesOfRow(node)) internalDofs.Add(dof);
                }
            }

            // The following are sorted in increasing order of boundary dof indices
            return (internalDofs.ToArray(), boundaryMultiplicities.Keys.ToArray(), boundaryMultiplicities.Values.ToArray(),
                boundaryConnectivities.Values.ToArray());
        }
    }
}
