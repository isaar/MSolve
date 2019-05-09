using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: perhaps the constrained dof orderings should not be optional. Instead there would be only 1 ISubdomainDofOrdering that 
//      handles both free and constrained dofs. Numbering constrained dofs is efficient and may not need to be repeated even if 
//      free dofs are renumbered (e.g. in XFEM).
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public abstract class SubdomainConstrainedDofOrderingBase : ISubdomainConstrainedDofOrdering
    {
        public SubdomainConstrainedDofOrderingBase(int numConstrainedDofs, DofTable subdomainConstrainedDofs)
        {
            this.NumConstrainedDofs = numConstrainedDofs;
            this.ConstrainedDofs = subdomainConstrainedDofs;
        }

        public DofTable ConstrainedDofs { get; }
        public int NumConstrainedDofs { get; }

        //TODO: if this class is merged with the ISubdomainFreeOrdering, this method should be accessed by the subdomain's 
        //      instance of ISubdomainFreeOrdering, instead of being static. It should  perform optimizations, such as using a 
        //      prescribed displacements vector and element constrained dof maps, similarly to how free displacements are 
        //      handled by ISubdomainFreeDofOrdering.
        public static void ApplyConstraintDisplacements(IElement element, double[] elementNodalDisplacements,
            Table<INode, IDofType, double> constraints)
        {
            int elementDofIdx = 0;
            IReadOnlyList<INode> nodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            IReadOnlyList<IReadOnlyList<IDofType>> dofs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            for (int i = 0; i < nodes.Count; ++i)
            {
                //bool isConstrainedNode = constraintsDictionary.TryGetValue(nodes[i].ID, 
                //    out Dictionary<DOFType, double> constrainedDOFs);
                bool isConstrainedNode = constraints.TryGetDataOfRow(nodes[i],
                    out IReadOnlyDictionary<IDofType, double> constrainedDOFs);
                if (isConstrainedNode)
                {
                    foreach (IDofType dofType in dofs[i])
                    {
                        bool isConstrainedDof = constrainedDOFs.TryGetValue(dofType, out double constraintDisplacement);
                        //if (isConstrainedNode && isConstrainedDof)
                        if (isConstrainedDof)
                        {
                            Debug.Assert(elementNodalDisplacements[elementDofIdx] == 0); // TODO: and why is this an assumption?
                            elementNodalDisplacements[elementDofIdx] = constraintDisplacement;
                        }
                        ++elementDofIdx;
                    }
                }
                else elementDofIdx += dofs[i].Count;
            }
        }

        protected (int numAllElementDofs, int[] elementDofIndices, int[] subdomainDofIndices)
            ProcessConstrainedDofsOfElement(IElement element)
        {
            IReadOnlyList<INode> elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            IReadOnlyList<IReadOnlyList<IDofType>> elementDofs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);

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

        public abstract (int[] elementDofIndices, int[] subdomainDofIndices) 
            MapConstrainedDofsElementToSubdomain(IElement element);
    }
}
