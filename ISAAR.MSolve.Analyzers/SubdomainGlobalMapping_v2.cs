using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    public class SubdomainGlobalMapping_v2 : ISubdomainGlobalMapping_v2
    {
        private readonly ISubdomain_v2 subdomain;

        public SubdomainGlobalMapping_v2(ISubdomain_v2 subdomain)
        {
            this.subdomain = subdomain;
        }

        public void SplitGlobalVectorToSubdomain(IVectorView globalVector, IVector subdomainVector)
        {
            subdomain.DofOrdering.ExtractVectorSubdomainFromGlobal(subdomain, globalVector, subdomainVector);
        }

        public void SubdomainToGlobalVector(IVectorView subdomainVector, IVector globalVector)
         => subdomain.DofOrdering.AddVectorSubdomainToGlobal(subdomain, subdomainVector, globalVector);

        public void SubdomainToGlobalVectorMeanValue(IVectorView subdomainVector, IVector globalVector)
        {
            throw new NotImplementedException();
        }
    }
}

