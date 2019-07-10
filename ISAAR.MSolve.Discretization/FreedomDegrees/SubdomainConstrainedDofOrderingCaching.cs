using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: This should be thread safe.
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class SubdomainConstrainedDofOrderingCaching : SubdomainConstrainedDofOrderingBase
    {
        private readonly Dictionary<IElement, (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices)> 
            elementDofsCache = new Dictionary<IElement, (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices)>();

        public SubdomainConstrainedDofOrderingCaching(int numConstrainedDofs, DofTable subdomainConstrainedDofs) :
            base(numConstrainedDofs, subdomainConstrainedDofs) { }

        public override (int[] elementDofIndices, int[] subdomainDofIndices) 
            MapConstrainedDofsElementToSubdomain(IElement element)
        {
            (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices) = GetElementData(element);
            return (elementDofIndices, subdomainDofIndices);
        }

        private (int numAllDofs, int[] elementDofIndices, int[] subdomainDofIndices) GetElementData(IElement element)
        {
            bool isStored = elementDofsCache.TryGetValue(element, out (int, int[], int[]) elementData);
            if (isStored) return elementData;
            else
            {
                elementData = base.ProcessConstrainedDofsOfElement(element);
                elementDofsCache.Add(element, elementData);
                return elementData;
            }
        }
    }
}
