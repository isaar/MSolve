using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    public class SubdomainGlobalMapping : ISubdomainGlobalMapping
    {
        private readonly Subdomain subdomain;
        public SubdomainGlobalMapping(Subdomain subdomain)
        {
            this.subdomain = subdomain;
        }

        public void SplitGlobalVectorToSubdomain(double[] vIn, double[] vOut)
        {
            this.subdomain.SplitGlobalVectorToSubdomain(vIn, vOut);
        }

        public void SubdomainToGlobalVector(double[] vIn, double[] vOut)
        {
            this.subdomain.SubdomainToGlobalVector(vIn, vOut);
        }

        public void SubdomainToGlobalVectorMeanValue(double[] vIn, double[] vOut)
        {
            throw new NotImplementedException();
        }
    }
}

