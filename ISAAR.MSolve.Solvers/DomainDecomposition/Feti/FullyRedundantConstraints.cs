using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    /// <summary>
    /// Applies fully redundant constraints
    /// </summary>
    internal class FullyRedundantConstraints : ICrosspointStrategy
    {
        public (Subdomain_v2[] subdomainsPlus, Subdomain_v2[] subdomainsMinus) FindSubdomainCombinations(int nodeMultiplicity,
            IEnumerable<Subdomain_v2> nodeSubdomains)
        {
            Debug.Assert(nodeMultiplicity > 1);
            int numNodeCombos = (nodeMultiplicity * (nodeMultiplicity - 1)) / 2; //TODO: not sure about this
            var subdomainsPlus = new Subdomain_v2[numNodeCombos];
            var subdomainsMinus = new Subdomain_v2[numNodeCombos];

            var processedSubdomains = new HashSet<Subdomain_v2>(nodeSubdomains);
            int comboCounter = 0;
            foreach (Subdomain_v2 subdomain1 in nodeSubdomains)
            {
                processedSubdomains.Remove(subdomain1);
                foreach (Subdomain_v2 subdomain2 in processedSubdomains)
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
