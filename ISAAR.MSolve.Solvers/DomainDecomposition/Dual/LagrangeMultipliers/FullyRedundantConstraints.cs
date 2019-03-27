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
        public (ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus) FindSubdomainCombinations(int nodeMultiplicity,
            IEnumerable<ISubdomain_v2> nodeSubdomains)
        {
            Debug.Assert(nodeMultiplicity > 1);
            int numNodeCombos = (nodeMultiplicity * (nodeMultiplicity - 1)) / 2; //TODO: not sure about this
            var subdomainsPlus = new ISubdomain_v2[numNodeCombos];
            var subdomainsMinus = new ISubdomain_v2[numNodeCombos];

            var processedSubdomains = new HashSet<ISubdomain_v2>(nodeSubdomains);
            int comboCounter = 0;
            foreach (ISubdomain_v2 subdomain1 in nodeSubdomains)
            {
                processedSubdomains.Remove(subdomain1);
                foreach (ISubdomain_v2 subdomain2 in processedSubdomains)
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
