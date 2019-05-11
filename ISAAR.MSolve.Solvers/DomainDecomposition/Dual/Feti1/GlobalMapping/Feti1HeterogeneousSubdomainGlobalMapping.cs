using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.GlobalMapping
{
    public class Feti1HeterogeneousSubdomainGlobalMapping : Feti1SubdomainGlobalMappingBase
    {
        private readonly HeterogeneousStiffnessDistribution distribution;

        public Feti1HeterogeneousSubdomainGlobalMapping(IStructuralModel model, Feti1DofSeparator dofSeparator,
            HeterogeneousStiffnessDistribution distribution) : base(model, dofSeparator)
        {
            this.distribution = distribution;
        }

        public override Dictionary<int, SparseVector> DistributeNodalLoads(
            IReadOnlyDictionary<int, ILinearSystem> linearSystems, Table<INode, IDofType, double> globalNodalLoads)
        {
            //TODO: This should be done using Dictionary<int, double[]> relativeBoundaryStiffnesses, instead of recreating that data.
            //TODO: Should I implement this as fb(s) = Lpb(s) * fb, Lpb(s) = Db(s)*Lb(s) * inv(Lb^T*Db*Lb)?
            //TODO: Internal loaded dofs should be handled differently as an optimization.

            var subdomainLoads = new Dictionary<int, SortedDictionary<int, double>>();
            foreach (var subdomainID in linearSystems.Keys) subdomainLoads[subdomainID] = new SortedDictionary<int, double>();

            foreach ((INode node, IDofType dofType, double loadAmount) in globalNodalLoads)
            {
                if (node.SubdomainsDictionary.Count == 1) // optimization for internal dof
                {
                    ISubdomain subdomain = node.SubdomainsDictionary.First().Value;
                    int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                    subdomainLoads[subdomain.ID][subdomainDofIdx] = loadAmount;
                }
                else // boundary dof: regularize with respect to the diagonal entries of the stiffness matrix at this dof
                {
                    BoundaryDofLumpedStiffness dofStiffness = distribution.BoundaryDofStiffnesses[node, dofType];
                    foreach (var idSubdomain in node.SubdomainsDictionary)
                    {
                        int id = idSubdomain.Key;
                        ISubdomain subdomain = idSubdomain.Value;
                        int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                        double relativeStiffness = dofStiffness.SubdomainStiffnesses[subdomain] / dofStiffness.TotalStiffness;
                        subdomainLoads[id][subdomainDofIdx] = loadAmount * relativeStiffness;
                    }
                }
            }

            return BuildForceVectors(linearSystems, subdomainLoads);
        }

        protected override double[] CalcBoundaryDofMultipliers(ISubdomain subdomain)
        {
            //TODO: Should this be cached? It stores the same info as HeterogeneousStiffnessDistribution.BoundaryDofStiffnesses.
            //      This format is more compact and has more efficient indexing when dealing with only 1 subdomain at a time, 
            //      but is difficult to index when accessing global boundary dofs. Is it possible to only use one of the two 
            //      formats? 
            //TODO: Could this be handled when extracting the lumped boundary stiffnesses? That way we can avoid searching
            //      the subdomain in BoundaryDofLumpedStiffness.SubdomainStiffnesses for each dof.

            int[] boundaryDofIndices = dofSeparator.BoundaryDofIndices[subdomain.ID];
            (INode, IDofType)[] boundaryDofs = dofSeparator.BoundaryDofs[subdomain.ID];
            int numBoundaryDofs = boundaryDofIndices.Length;
            var relativeStiffnesses = new double[numBoundaryDofs];
            for (int i = 0; i < boundaryDofIndices.Length; ++i)
            {
                (INode node, IDofType dofType) = boundaryDofs[i];
                BoundaryDofLumpedStiffness dofStiffness = distribution.BoundaryDofStiffnesses[node, dofType];
                double relativeStiffness = dofStiffness.SubdomainStiffnesses[subdomain] / dofStiffness.TotalStiffness;
                relativeStiffnesses[i] = relativeStiffness;
            }
            return relativeStiffnesses;
        }
    }
}
