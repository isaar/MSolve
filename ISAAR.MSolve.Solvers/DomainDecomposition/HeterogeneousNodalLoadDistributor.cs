using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: Should I implement this as fb(s) = Lpb(s) * fb, Lpb(s) = Db(s)*Lb(s) * inv(Lb^T*Db*Lb)?
//TODO: Internal loaded dofs should be handled differently as an optimization.
namespace ISAAR.MSolve.Solvers.DomainDecomposition
{
    public class HeterogeneousNodalLoadDistributor : INodalLoadDistributor
    {
        public Dictionary<int, SparseVector> DistributeNodalLoads(IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems,
            Table<INode, DOFType, double> globalNodalLoads)
        {
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

        private Dictionary<int, SparseVector> BuildForceVectors(IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems,
            Dictionary<int, SortedDictionary<int, double>> subdomainLoads)
        {
            var result = new Dictionary<int, SparseVector>();
            foreach (var idLoads in subdomainLoads)
            {
                int numSubdomainDofs = linearSystems[idLoads.Key].Subdomain.FreeDofOrdering.NumFreeDofs;
                result[idLoads.Key] = SparseVector.CreateFromDictionary(numSubdomainDofs, idLoads.Value);
            }
            return result;
        }
    }
}
