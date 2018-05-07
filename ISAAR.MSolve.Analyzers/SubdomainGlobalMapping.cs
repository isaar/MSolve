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
            foreach (int nodeID in this.subdomain.GlobalNodalDOFsDictionary.Keys)
            {
                Dictionary<DOFType, int> dofTypes = this.subdomain.NodalDOFsDictionary[nodeID];
                foreach (DOFType dofType in dofTypes.Keys)
                {
                    int localDOF = this.subdomain.NodalDOFsDictionary[nodeID][dofType];
                    int globalDOF = this.subdomain.GlobalNodalDOFsDictionary[nodeID][dofType];
                    if (localDOF > -1 && globalDOF > -1) vOut[localDOF] = vIn[globalDOF];
                }
            }
        }

        public void SubdomainToGlobalVector(double[] vIn, double[] vOut)
        {
            foreach (int nodeID in this.subdomain.GlobalNodalDOFsDictionary.Keys)
            {
                Dictionary<DOFType, int> dofTypes = this.subdomain.NodalDOFsDictionary[nodeID];
                foreach (DOFType dofType in dofTypes.Keys)
                {
                    int localDOF = this.subdomain.NodalDOFsDictionary[nodeID][dofType];
                    int globalDOF = this.subdomain.GlobalNodalDOFsDictionary[nodeID][dofType];
                    if (localDOF > -1 && globalDOF > -1) vOut[globalDOF] += vIn[localDOF];
                }
            }
        }

        public void SubdomainToGlobalVectorMeanValue(double[] vIn, double[] vOut)
        {
            throw new NotImplementedException();
        }
    }
}
