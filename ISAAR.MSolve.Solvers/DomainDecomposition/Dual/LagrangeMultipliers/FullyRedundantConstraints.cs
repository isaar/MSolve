using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    /// <summary>
    /// Applies fully redundant constraints
    /// </summary>
    internal class FullyRedundantConstraints : ICrosspointStrategy
    {
        public (ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus) FindSubdomainCombinations(ISubdomain[] nodeSubdomains)
        {
            int nodeMultiplicity = nodeSubdomains.Length;
            Debug.Assert(nodeMultiplicity > 1);
            int numNodeCombos = (nodeMultiplicity * (nodeMultiplicity - 1)) / 2; //TODO: not sure about this
            var subdomainsPlus = new ISubdomain[numNodeCombos];
            var subdomainsMinus = new ISubdomain[numNodeCombos];

            var processedSubdomains = new HashSet<ISubdomain>(nodeSubdomains);
            int comboCounter = 0;
            foreach (ISubdomain subdomain1 in nodeSubdomains)
            {
                processedSubdomains.Remove(subdomain1);
                foreach (ISubdomain subdomain2 in processedSubdomains)
                {
                    subdomainsPlus[comboCounter] = subdomain1;
                    subdomainsMinus[comboCounter] = subdomain2;
                    ++comboCounter;
                }
            }
            return (subdomainsPlus, subdomainsMinus);
        }

    }
}
