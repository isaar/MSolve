using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    public class SubdomainGlobalMapping_v2 : ISubdomainGlobalMapping_v2
    {
        private readonly Subdomain_v2 subdomain;

        public SubdomainGlobalMapping_v2(Subdomain_v2 subdomain)
        {
            this.subdomain = subdomain;
        }

        public void SplitGlobalVectorToSubdomain(IVectorView globalVector, IVector subdomainVector)
        {
            //subdomain.SplitGlobalVectorToSubdomain_v2(globalVector, subdomainVector);
            subdomain.SplitGlobalVectorToSubdomain_v3(globalVector, subdomainVector);
        }

        public void SubdomainToGlobalVector(IVectorView subdomainVector, IVector globalVector)
        {
            //this.subdomain.SubdomainToGlobalVector_v2(subdomainVector, globalVector);
            subdomain.SubdomainToGlobalVector_v3(subdomainVector, globalVector);
        }

        public void SubdomainToGlobalVectorMeanValue(IVectorView subdomainVector, IVector globalVector)
        {
            throw new NotImplementedException();
        }
    }
}

