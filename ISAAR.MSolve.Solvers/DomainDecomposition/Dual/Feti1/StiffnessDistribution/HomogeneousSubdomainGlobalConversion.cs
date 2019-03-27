using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution
{
    internal class HomogeneousSubdomainGlobalConversion : Feti1SubdomainGlobalConversionBase
    {
        internal HomogeneousSubdomainGlobalConversion(IStructuralModel_v2 model, Feti1DofSeparator dofSeparator) :
            base(model, dofSeparator)
        {
        }

        public override Dictionary<int, SparseVector> DistributeNodalLoads(
            IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems, Table<INode, DOFType, double> globalNodalLoads)
        {
            //TODO: Should I implemented this as fb(s) = Lpb(s) * fb, Lpb(s) = Lb(s) * inv(Mb)?
            //TODO: Internal loaded dofs should be handled differently as an optimization.

            var subdomainLoads = new Dictionary<int, SortedDictionary<int, double>>();
            foreach (var subdomainID in linearSystems.Keys) subdomainLoads[subdomainID] = new SortedDictionary<int, double>();

            foreach ((INode node, DOFType dofType, double amount) in globalNodalLoads)
            {
                double amountPerSubdomain = amount / node.SubdomainsDictionary.Count;
                foreach (var idSubdomain in node.SubdomainsDictionary)
                {
                    int subdomainDofIdx = idSubdomain.Value.FreeDofOrdering.FreeDofs[node, dofType];
                    subdomainLoads[idSubdomain.Key][subdomainDofIdx] = amountPerSubdomain;
                }
            }

            return BuildForceVectors(linearSystems, subdomainLoads);
        }

        protected override double[] CalcBoundaryDofMultipliers(ISubdomain_v2 subdomain)
        {
            int[] multiplicites = dofSeparator.BoundaryDofMultiplicities[subdomain.ID];
            var inverseMultiplicites = new double[multiplicites.Length];
            for (int i = 0; i < multiplicites.Length; ++i) inverseMultiplicites[i] = 1.0 / multiplicites[i];
            return inverseMultiplicites;
        }
    }
}
