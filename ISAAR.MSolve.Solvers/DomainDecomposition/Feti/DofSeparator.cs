using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Needs a proper name. This probably cannot be incorporated in the ISubdomainDofOrdering, as the intent is different and
//      depending on the DD method the dof categories may be different (e.g. FETI-1: internal/boundary, 
//      FETI-DP: corner/boundary/remainder)
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public class DofSeparator
    {
        internal Dictionary<int, int[]> BoundaryDofs { get; private set; }
        internal Dictionary<int, int[]> BoundaryDofsMultiplicity { get; private set; }
        internal Dictionary<INode, DOFType[]> GlobalBoundaryDofs { get; private set; }
        internal Dictionary<int, int[]> InternalDofs { get; private set; }

        internal void SeparateBoundaryInternalDofs(IStructuralModel_v2 model)
        {
            GatherGlobalBoundaryDofs(model);
            SeparateBoundaryInternalDofsOfSubdomains(model);
        }

        private void SeparateBoundaryInternalDofsOfSubdomains(IStructuralModel_v2 model)
        {
            BoundaryDofs = new Dictionary<int, int[]>();
            BoundaryDofsMultiplicity = new Dictionary<int, int[]>();
            InternalDofs = new Dictionary<int, int[]>();
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                var boundaryDofsOfSubdomain = new SortedDictionary<int, int>(); // key = dofIdx, value = multiplicity
                var internalDofsOfSubdomain = new SortedSet<int>();
                foreach (INode node in subdomain.Nodes)
                {
                    IEnumerable<int> nodalDofs = subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node);
                    int nodeMultiplicity = node.SubdomainsDictionary.Count;
                    if (nodeMultiplicity > 1) // boundary node
                    {
                        foreach (int dof in nodalDofs) boundaryDofsOfSubdomain.Add(dof, nodeMultiplicity);
                    }
                    else
                    {
                        foreach (int dof in nodalDofs) internalDofsOfSubdomain.Add(dof);
                    }
                }
                BoundaryDofs.Add(subdomain.ID, boundaryDofsOfSubdomain.Keys.ToArray());
                BoundaryDofsMultiplicity.Add(subdomain.ID, boundaryDofsOfSubdomain.Values.ToArray()); // sorted the same as Keys
                InternalDofs.Add(subdomain.ID, internalDofsOfSubdomain.ToArray());
            }
        }

        //TODO: This has common code with the LagrangeMultipliersEnumerator.DefineBoundaryDofConstraints(). Create the object 
        //      from here, store it somewhere (could be here) and pass the global boundary dofs to that method.
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
