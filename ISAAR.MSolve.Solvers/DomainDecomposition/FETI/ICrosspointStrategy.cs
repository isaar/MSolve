using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: extend this to methods other than FEM.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.FETI
{
    /// <summary>
    /// Finds the combinations of subdomains that contain same boundary node.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal interface ICrosspointStrategy
    {
        //TODO: Perhaps a list of pairs is better than a pair of lists.
        (Subdomain_v2[] subdomainsPlus, Subdomain_v2[] subdomainsMinus) FindSubdomainCombinations(int nodeMultiplicity,
            IEnumerable<Subdomain_v2> nodeSubdomains);
    }
}
