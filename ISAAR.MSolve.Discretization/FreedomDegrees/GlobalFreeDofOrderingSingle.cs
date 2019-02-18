using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class GlobalFreeDofOrderingSingle: IGlobalFreeDofOrdering
    {
        private readonly ISubdomain_v2 subdomain;
        private readonly int[] subdomainToGlobalDofMap;

        public GlobalFreeDofOrderingSingle(ISubdomain_v2 subdomain, ISubdomainFreeDofOrdering subdomainDofOrdering)
        {
            this.subdomain = subdomain;

            this.NumGlobalFreeDofs = subdomainDofOrdering.NumFreeDofs;
            this.GlobalFreeDofs = subdomainDofOrdering.FreeDofs;
            this.SubdomainDofOrderings = 
                new Dictionary<ISubdomain_v2, ISubdomainFreeDofOrdering> { { subdomain, subdomainDofOrdering } };

            // A for-loop is faster than LINQ
            this.subdomainToGlobalDofMap = new int[NumGlobalFreeDofs];
            for (int i = 0; i < NumGlobalFreeDofs; ++i) this.subdomainToGlobalDofMap[i] = i;
        }

        public DofTable GlobalFreeDofs { get; }
        public int NumGlobalFreeDofs { get; }
        public IReadOnlyDictionary<ISubdomain_v2, ISubdomainFreeDofOrdering> SubdomainDofOrderings { get; }

        //TODO: Usually the only vectors added are the ones from subdomains, therefore in this case we can just copy
        //      subdomainVector, instead of adding it.
        public void AddVectorSubdomainToGlobal(ISubdomain_v2 subdomain, IVectorView subdomainVector, IVector globalVector)
        {
            Debug.Assert(subdomain == this.subdomain);
            globalVector.AddSubvectorIntoThis(0, subdomainVector, 0, subdomainVector.Length);
        }

        public void AddVectorSubdomainToGlobalMeanValue(ISubdomain_v2 subdomain, IVectorView subdomainVector, 
            IVector globalVector) => throw new NotImplementedException();

        public void ExtractVectorSubdomainFromGlobal(ISubdomain_v2 subdomain, IVectorView globalVector, IVector subdomainVector)
        {
            Debug.Assert(subdomain == this.subdomain);
            subdomainVector.CopyFrom(globalVector);
        }

        //TODO: Using this map is a waste of computational power: map[i] = i
        public int[] MapFreeDofsSubdomainToGlobal(ISubdomain_v2 subdomain)
        {
            Debug.Assert(subdomain == this.subdomain);
            return subdomainToGlobalDofMap;
        }
    }
}
