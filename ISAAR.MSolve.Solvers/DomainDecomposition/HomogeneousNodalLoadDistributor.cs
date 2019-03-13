using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: Should I implement this as fb(s) = Lpb(s) * fb, Lpb(s) = Lb(s) * inv(Mb)?
//TODO: Internal loaded dofs should be handled differently as an optimization.
namespace ISAAR.MSolve.Solvers.DomainDecomposition
{
    public class HomogeneousNodalLoadDistributor : INodalLoadDistributor
    {
        public Dictionary<int, SparseVector> DistributeNodalLoads(IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems,
            Table<INode, DOFType, double> globalNodalLoads)
        {
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
