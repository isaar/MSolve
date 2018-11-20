using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: perhaps I should store the subdomain orderings as well.
//TODO: implement a global ordering optimized for the case where there is only 1 subdomain.
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class GlobalFreeDofOrderingGeneral: IGlobalFreeDofOrdering
    {
        private readonly Dictionary<ISubdomain_v2, int[]> subdomainToGlobalDofMaps;

        public GlobalFreeDofOrderingGeneral(int numGlobalFreeDofs, DofTable globalFreeDofs, 
            Dictionary<ISubdomain_v2, ISubdomainFreeDofOrdering> subdomainDofOrderings)
        {
            this.NumGlobalFreeDofs = numGlobalFreeDofs;
            this.GlobalFreeDofs = globalFreeDofs;
            this.SubdomainDofOrderings = subdomainDofOrderings;

            subdomainToGlobalDofMaps = new Dictionary<ISubdomain_v2, int[]>(subdomainDofOrderings.Count);
            foreach (var subdomainOrderingPair in subdomainDofOrderings)
            {
                var subdomainToGlobalDofs = new int[subdomainOrderingPair.Value.NumFreeDofs];
                foreach ((INode node, DOFType dofType, int subdomainDofIdx) in subdomainOrderingPair.Value.FreeDofs)
                {
                    subdomainToGlobalDofs[subdomainDofIdx] = globalFreeDofs[node, dofType];
                }
                subdomainToGlobalDofMaps.Add(subdomainOrderingPair.Key, subdomainToGlobalDofs);
            }
        }

        public DofTable GlobalFreeDofs { get; }
        public int NumGlobalFreeDofs { get; }
        public IReadOnlyDictionary<ISubdomain_v2, ISubdomainFreeDofOrdering> SubdomainDofOrderings { get; }

        public void AddVectorSubdomainToGlobal(ISubdomain_v2 subdomain, IVectorView subdomainVector, IVector globalVector)
        {
            ISubdomainFreeDofOrdering subdomainOrdering = SubdomainDofOrderings[subdomain];
            int[] subdomainToGlobalDofs = subdomainToGlobalDofMaps[subdomain];

            for (int i = 0; i < subdomainOrdering.NumFreeDofs; ++i)
            {
                int globalIdx = subdomainToGlobalDofs[i];
                globalVector.Set(globalIdx, subdomainVector[i] + globalVector[globalIdx]);
                //TODO: add a Vector.SetSubvector and Vector.AddSubvector for incontiguous entries
            }
        }

        public void AddVectorSubdomainToGlobalMeanValue(ISubdomain_v2 subdomain, IVectorView subdomainVector, 
            IVector globalVector) => throw new NotImplementedException();

        public void ExtractVectorSubdomainFromGlobal(ISubdomain_v2 subdomain, IVectorView globalVector, IVector subdomainVector)
        {
            ISubdomainFreeDofOrdering subdomainOrdering = SubdomainDofOrderings[subdomain];
            int[] subdomainToGlobalDofs = subdomainToGlobalDofMaps[subdomain];

            for (int i = 0; i < subdomainOrdering.NumFreeDofs; ++i)
            {
                //TODO: add a Vector.SetSubvector and Vector.AddSubvector for incontiguous entries
                subdomainVector.Set(i, globalVector[subdomainToGlobalDofs[i]]);
            }
        }

        //TODO: the returned array should be readonly.
        public int[] MapFreeDofsSubdomainToGlobal(ISubdomain_v2 subdomain) => subdomainToGlobalDofMaps[subdomain];
    }
}
