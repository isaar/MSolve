using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: This should be a static class that will be called by the DofSeparator of each FETI method to create the common data structures.
//TODO: Should this be abstract or not suffixed with "Base".
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual
{
    public abstract class DofSeparatorBase : IDofSeparator
    {
        public abstract Dictionary<int, int[]> BoundaryDofIndices { get; protected set; }
        public abstract Dictionary<int, (INode node, IDofType dofType)[]> BoundaryDofs { get; protected set; }
        public Dictionary<INode, IDofType[]> GlobalBoundaryDofs { get; private set; }
        public abstract Dictionary<int, int[]> InternalDofIndices { get; protected set; }

        protected void DefineGlobalBoundaryDofs(IEnumerable<INode> allNodes, IGlobalFreeDofOrdering globalDofOrdering)
        {
            GlobalBoundaryDofs = new Dictionary<INode, IDofType[]>();

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
                    else GlobalBoundaryDofs[node] = dofsOfNode;
                }
            }
        }

        protected static (int[] internalDofIndices, int[] boundaryDofIndices, (INode node, IDofType dofType)[] boundaryDofs) 
            SeparateBoundaryInternalDofs(IEnumerable<INode> nodes, DofTable freeDofs)
        {
            var boundaryDofs = new SortedDictionary<int, (INode node, IDofType dofType)>(); // key = dofIdx, value = (node, dofType)
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
                            boundaryDofs.Add(dofTypeIdxPair.Value, (node, dofTypeIdxPair.Key));
                        }
                    }
                }
                else // internal nodes
                {
                    foreach (int dof in freeDofs.GetValuesOfRow(node)) internalDofs.Add(dof);
                }
            }

            // The following are sorted in increasing order of boundary dof indices
            return (internalDofs.ToArray(), boundaryDofs.Keys.ToArray(), boundaryDofs.Values.ToArray());
        }
    }
}
