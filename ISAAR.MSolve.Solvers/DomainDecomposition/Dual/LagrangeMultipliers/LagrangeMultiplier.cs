using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: perhaps I should also store the dof indices in each subdomain. In that case, creating the boolean matrices can be 
//      decoupled from the lagrange enumeration.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    internal class LagrangeMultiplier
    {
        internal LagrangeMultiplier(INode node, DOFType dof, ISubdomain_v2 subdomainPlus, ISubdomain_v2 subdomainMinus)
        {
            this.Node = node;
            this.DofType = dof;
            this.SubdomainPlus = subdomainPlus;
            this.SubdomainMinus = subdomainMinus;
        }

        internal DOFType DofType { get; }
        internal INode Node { get; }
        internal ISubdomain_v2 SubdomainMinus { get; }
        internal ISubdomain_v2 SubdomainPlus { get; }
    }
}
