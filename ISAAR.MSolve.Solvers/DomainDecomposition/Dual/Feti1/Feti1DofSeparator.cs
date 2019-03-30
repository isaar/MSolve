using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Needs a proper name. This probably cannot be incorporated in the ISubdomainDofOrdering, as the intent is different and
//      depending on the DD method the dof categories may be different (e.g. FETI-1: internal/boundary, 
//      FETI-DP: corner/boundary/remainder)
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1
{
    public class Feti1DofSeparator : IDofSeparator
    {
        //TODO: BoundaryDofConnectivities is useful only for heterogeneous problems and its creation might be costly.
        //      Needs benchmarking and possibly a way to avoid it in homogeneous problems
        public Dictionary<INode, DOFType[]> GlobalBoundaryDofs { get; private set; }
        internal Dictionary<int, (INode node, DOFType dofType)[]> BoundaryDofConnectivities { get; private set; }
        internal Dictionary<int, int[]> BoundaryDofIndices { get; private set; }
        internal Dictionary<int, int[]> BoundaryDofMultiplicities { get; private set; }
        internal Dictionary<int, int[]> InternalDofIndices { get; private set; }

        internal void SeparateBoundaryInternalDofs(IStructuralModel_v2 model)
        {
            GatherGlobalBoundaryDofs(model);
            SeparateBoundaryInternalDofsOfSubdomains(model);
        }

        private void SeparateBoundaryInternalDofsOfSubdomains(IStructuralModel_v2 model)
        {
            BoundaryDofIndices = new Dictionary<int, int[]>();
            BoundaryDofMultiplicities = new Dictionary<int, int[]>();
            BoundaryDofConnectivities = new Dictionary<int, (INode node, DOFType dofType)[]>();
            InternalDofIndices = new Dictionary<int, int[]>();
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                var boundaryMultiplicitiesOfSubdomain = new SortedDictionary<int, int>(); // key = dofIdx, value = multiplicity
                var boundaryConnectivitiesOfSubdomain = new SortedDictionary<int, (INode node, DOFType dofType)>(); 
                var internalDofsOfSubdomain = new SortedSet<int>();
                foreach (INode node in subdomain.Nodes)
                {
                    IEnumerable<int> nodalDofs = subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node);
                    int nodeMultiplicity = node.SubdomainsDictionary.Count;
                    if (nodeMultiplicity > 1) // boundary node
                    {
                        bool isNodeFree = subdomain.FreeDofOrdering.FreeDofs.TryGetDataOfRow(node, 
                            out IReadOnlyDictionary<DOFType, int> dofTypesIndices); // This avoids embedded and constrained dofs
                        if (isNodeFree)
                        {
                            foreach (var dofTypeIdxPair in dofTypesIndices)
                            {
                                boundaryMultiplicitiesOfSubdomain.Add(dofTypeIdxPair.Value, nodeMultiplicity);
                                boundaryConnectivitiesOfSubdomain.Add(dofTypeIdxPair.Value, (node, dofTypeIdxPair.Key));
                            }
                        }
                    }
                    else // internal nodes
                    {
                        foreach (int dof in subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node))
                        {
                            internalDofsOfSubdomain.Add(dof);
                        }
                    }
                }

                // The following are sorted in increasing order of boundary dof indices
                BoundaryDofIndices.Add(subdomain.ID, boundaryMultiplicitiesOfSubdomain.Keys.ToArray());
                BoundaryDofMultiplicities.Add(subdomain.ID, boundaryMultiplicitiesOfSubdomain.Values.ToArray());
                BoundaryDofConnectivities.Add(subdomain.ID, boundaryConnectivitiesOfSubdomain.Values.ToArray());

                InternalDofIndices.Add(subdomain.ID, internalDofsOfSubdomain.ToArray());
            }
        }

        private void GatherGlobalBoundaryDofs(IStructuralModel_v2 model)
        {
            GlobalBoundaryDofs = new Dictionary<INode, DOFType[]>();

            //TODO: model.Nodes probably doesn't work if there are embedded nodes. It is time to isolate the embedded nodes. Or I could use the GlobalDofOrdering.
            foreach (INode node in model.Nodes) 
            {
                int nodeMultiplicity = node.SubdomainsDictionary.Count;
                if (nodeMultiplicity > 1)
                {
                    // Access the free dofs only. Does this also filter out embedded dofs?
                    DOFType[] dofsOfNode = model.GlobalDofOrdering.GlobalFreeDofs.GetColumnsOfRow(node).ToArray(); //TODO: interacting with the table can be optimized.

                    // If all dofs of this node are constrained, then it is not considered boundary.
                    if (dofsOfNode.Length == 0) continue;
                    else GlobalBoundaryDofs[node] = dofsOfNode;
                }
            }
        }
    }
}
