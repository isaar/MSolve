using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution
{
    internal class HeterogeneousSubdomainGlobalConversion : Feti1SubdomainGlobalConversionBase
    {
        private readonly Dictionary<int, double[]> relativeBoundaryStiffnesses;

        internal HeterogeneousSubdomainGlobalConversion(IStructuralModel_v2 model, Feti1DofSeparator dofSeparator,
            Dictionary<int, double[]> relativeBoundaryStiffnesses) : base(model, dofSeparator)
        {
            this.relativeBoundaryStiffnesses = relativeBoundaryStiffnesses;
        }

        public override Dictionary<int, SparseVector> DistributeNodalLoads(
            IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems, Table<INode, DOFType, double> globalNodalLoads)
        {
            //TODO: This should be done using Dictionary<int, double[]> relativeBoundaryStiffnesses, instead of recreating that data.
            //TODO: Should I implement this as fb(s) = Lpb(s) * fb, Lpb(s) = Db(s)*Lb(s) * inv(Lb^T*Db*Lb)?
            //TODO: Internal loaded dofs should be handled differently as an optimization.

            var subdomainLoads = new Dictionary<int, SortedDictionary<int, double>>();
            foreach (var subdomainID in linearSystems.Keys) subdomainLoads[subdomainID] = new SortedDictionary<int, double>();

            foreach ((INode node, DOFType dofType, double amount) in globalNodalLoads)
            {
                if (node.SubdomainsDictionary.Count == 1) // optimization for internal dof
                {
                    ISubdomain_v2 subdomain = node.SubdomainsDictionary.First().Value;
                    int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                    subdomainLoads[subdomain.ID][subdomainDofIdx] = amount;
                }
                else // boundary dof: regularize with respect to the diagonal entries of the stiffness matrix at this dof
                {
                    ISubdomain_v2[] dofSubdomains = node.SubdomainsDictionary.Values.ToArray();
                    var dofIndices = new int[dofSubdomains.Length];
                    var dofStiffnesses = new double[dofSubdomains.Length];
                    for (int s = 0; s < dofSubdomains.Length; ++s)
                    {
                        ISubdomain_v2 subdomain = dofSubdomains[s];
                        dofIndices[s] = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                        dofStiffnesses[s] = linearSystems[subdomain.ID].Matrix[dofIndices[s], dofIndices[s]]; //TODO: accessing diagonal entries should be optimized (and batched).
                        Debug.Assert(dofStiffnesses[s] > 0);
                    }

                    double loadOverStiffnessSum = amount / dofStiffnesses.Sum();
                    for (int s = 0; s < dofSubdomains.Length; ++s)
                    {
                        subdomainLoads[dofSubdomains[s].ID][dofIndices[s]] = loadOverStiffnessSum * dofStiffnesses[s];
                    }
                }
            }

            return BuildForceVectors(linearSystems, subdomainLoads);
        }

        protected override double[] CalcBoundaryDofMultipliers(ISubdomain_v2 subdomain)
            => relativeBoundaryStiffnesses[subdomain.ID];
    }
}
