using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Ideally all the relevant methods should return Vector (or at least Vector for element level) 
//TODO: The map element to subdomain should be a BiList, so that it can be read in both directions and passed to native code 
//      (e.g. CUDA). It should also be cached for each element, unless there isn't enough memory or if the dof ordering is 
//      updated frequently. It can also be used for mapping vectors (e.g. ExtractVectorElementFromSubdomain), not only the
//      stiffness matrix.
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public interface ISubdomainFreeDofOrdering
    {
        DofTable FreeDofs { get; }

        int NumFreeDofs { get; }

        //TODO: What should it contain for constrained dofs?
        void AddVectorElementToSubdomain(IElement element, IVectorView elementVector, IVector subdomainVector);

        int CountElementDofs(IElement element);

        //TODO: What should it contain for constrained dofs?
        //TODOMaria: here is where the element displacements are assigned to zero if they are constrained
        void ExtractVectorElementFromSubdomain(IElement element, IVectorView subdomainVector, IVector elementVector); 

        (int[] elementDofIndices, int[] subdomainDofIndices) MapFreeDofsElementToSubdomain(IElement element);

        //TODO: perhaps the subdomain should be passed in the constructor.
        void Reorder(IReorderingAlgorithm reorderingAlgorithm, ISubdomain_v2 subdomain);

        void ReorderNodeMajor(IReadOnlyList<INode> sortedNodes);
    }
}
