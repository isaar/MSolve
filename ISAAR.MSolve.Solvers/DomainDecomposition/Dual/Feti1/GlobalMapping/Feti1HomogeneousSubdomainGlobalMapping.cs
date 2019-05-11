using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.GlobalMapping
{
    public class Feti1HomogeneousSubdomainGlobalMapping : Feti1SubdomainGlobalMappingBase
    {
        private readonly HomogeneousStiffnessDistribution distribution;

        public Feti1HomogeneousSubdomainGlobalMapping(IStructuralModel model, Feti1DofSeparator dofSeparator,
            HomogeneousStiffnessDistribution distribution) : base(model, dofSeparator)
        {
            this.distribution = distribution;
        }

        public override Dictionary<int, SparseVector> DistributeNodalLoads(
            IReadOnlyDictionary<int, ILinearSystem> linearSystems, Table<INode, IDofType, double> globalNodalLoads)
        {
            //TODO: Should I implemented this as fb(s) = Lpb(s) * fb, Lpb(s) = Lb(s) * inv(Mb)?
            //TODO: Internal loaded dofs should be handled differently as an optimization.

            var subdomainLoads = new Dictionary<int, SortedDictionary<int, double>>();
            foreach (var subdomainID in linearSystems.Keys) subdomainLoads[subdomainID] = new SortedDictionary<int, double>();

            foreach ((INode node, IDofType dofType, double loadAmount) in globalNodalLoads)
            {
                int multiplicity = node.SubdomainsDictionary.Count;
                if (multiplicity == 1) // optimization for internal dof
                {
                    ISubdomain subdomain = node.SubdomainsDictionary.First().Value;
                    int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                    subdomainLoads[subdomain.ID][subdomainDofIdx] = loadAmount;
                }
                else
                {
                    double amountPerSubdomain = loadAmount / multiplicity;
                    foreach (var idSubdomain in node.SubdomainsDictionary)
                    {
                        int subdomainDofIdx = idSubdomain.Value.FreeDofOrdering.FreeDofs[node, dofType];
                        subdomainLoads[idSubdomain.Key][subdomainDofIdx] = amountPerSubdomain;
                    }
                }
            }

            return BuildForceVectors(linearSystems, subdomainLoads);
        }

        protected override double[] CalcBoundaryDofMultipliers(ISubdomain subdomain)
        {
            int[] multiplicites = distribution.BoundaryDofMultiplicities[subdomain.ID];
            var inverseMultiplicites = new double[multiplicites.Length];
            for (int i = 0; i < multiplicites.Length; ++i) inverseMultiplicites[i] = 1.0 / multiplicites[i];
            return inverseMultiplicites;
        }
    }
}
