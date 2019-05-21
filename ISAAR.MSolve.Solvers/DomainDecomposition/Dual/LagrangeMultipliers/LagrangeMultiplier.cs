using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: perhaps I should also store the dof indices in each subdomain. In that case, creating the boolean matrices can be 
//      decoupled from the lagrange enumeration.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    public class LagrangeMultiplier
    {
        internal LagrangeMultiplier(INode node, IDofType dof, ISubdomain subdomainPlus, ISubdomain subdomainMinus)
        {
            this.Node = node;
            this.DofType = dof;
            this.SubdomainPlus = subdomainPlus;
            this.SubdomainMinus = subdomainMinus;
        }

        internal IDofType DofType { get; }
        internal INode Node { get; }
        internal ISubdomain SubdomainMinus { get; }
        internal ISubdomain SubdomainPlus { get; }
    }
}
