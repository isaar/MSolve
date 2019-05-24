using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    public class FetiDPSubdomainGlobalMapping
    {
        private readonly IStiffnessDistribution distribution;
        private readonly FetiDPDofSeparator dofSeparator;
        private readonly IStructuralModel model;

        public FetiDPSubdomainGlobalMapping(IStructuralModel model, FetiDPDofSeparator dofSeparator,
            IStiffnessDistribution distribution)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.distribution = distribution;
        }

        public double CalculateGlobalForcesNorm(Dictionary<int, IVectorView> subdomainForces)
        {
            //TODO: This can be optimized: calculate the dot product f*f for the internal dofs of each subdomain separately,
            //      only assemble global vector for the boundary dofs, find its dot product with itself, add the contributions
            //      for the internal dofs and finally apply SQRT(). This would greatly reduce the communication requirements.
            //TODO: this should be used for non linear analyzers as well (instead of building the global RHS)
            //TODO: Is this correct? For the residual, it would be wrong to find f-K*u for each subdomain and then call this.

            return GatherGlobalForces(subdomainForces).Norm2();
        }

        public Dictionary<int, SparseVector> DistributeNodalLoads(Dictionary<int, ISubdomain> subdomains,
            Table<INode, IDofType, double> globalNodalLoads)
        {
            //TODO: Should I implement this as fb(s) = Lpb(s) * fb, with a) Lpb(s) = Lb(s) * inv(Mb) for homogeneous and 
            //      b) Lpb(s) = Db(s)*Lb(s) * inv(Lb^T*Db*Lb) for heterogeneous?

            var subdomainLoads = new Dictionary<int, SortedDictionary<int, double>>();
            foreach (var subdomainID in subdomains.Keys) subdomainLoads[subdomainID] = new SortedDictionary<int, double>();

            foreach ((INode node, IDofType dofType, double loadAmount) in globalNodalLoads)
            {
                bool isCornerDof = dofSeparator.GlobalCornerDofOrdering.Contains(node, dofType);
                if (isCornerDof) 
                {
                    // Loads at corner dofs will be distributed equally. It shouldn't matter how I distribute these, since I 
                    // will only sum them together again during the static condensation of remainder dofs phase.
                    //TODO: is that correct?
                    double loadPerSubdomain = loadAmount / node.SubdomainsDictionary.Count;
                    foreach (var idSubdomain in node.SubdomainsDictionary)
                    {
                        int id = idSubdomain.Key;
                        ISubdomain subdomain = idSubdomain.Value;
                        int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                        subdomainLoads[id][subdomainDofIdx] = loadPerSubdomain;
                    }
                }
                else
                {
                    if (node.SubdomainsDictionary.Count == 1) // optimization for internal dof
                    {
                        ISubdomain subdomain = node.SubdomainsDictionary.First().Value;
                        int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                        subdomainLoads[subdomain.ID][subdomainDofIdx] = loadAmount;
                    }
                    else // boundary dof: regularize with respect to the diagonal entries of the stiffness matrix at this dof
                    {
                        Dictionary<int, double> boundaryDofCoeffs = distribution.CalcBoundaryDofCoefficients(node, dofType);
                        foreach (var idSubdomain in node.SubdomainsDictionary)
                        {
                            int id = idSubdomain.Key;
                            ISubdomain subdomain = idSubdomain.Value;
                            int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                            subdomainLoads[id][subdomainDofIdx] = loadAmount * boundaryDofCoeffs[id];
                        }
                    }
                }
            }

            var subdomainVectors = new Dictionary<int, SparseVector>();
            foreach (var idSubdomains in subdomains)
            {
                int id = idSubdomains.Key;
                int numSubdomainDofs = idSubdomains.Value.FreeDofOrdering.NumFreeDofs;
                subdomainVectors[id] = SparseVector.CreateFromDictionary(numSubdomainDofs, subdomainLoads[id]);
            }
            return subdomainVectors;
        }

        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainRemainderDisplacements, 
            IVectorView globalCornerDisplacements)
        {
            var globalDisplacements = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);

            // Remainder dofs
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                int[] freeToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                int[] remainderToFreeDofs = dofSeparator.RemainderDofIndices[id];
                IVectorView remainderDisplacements = subdomainRemainderDisplacements[id]; //TODO: benchmark the performance if this was concrete Vector

                // Internal dofs are copied without averaging.
                foreach (int remainderDofIdx in dofSeparator.InternalDofIndices[id])
                {
                    int globalDofIdx = freeToGlobalDofs[remainderToFreeDofs[remainderDofIdx]];
                    globalDisplacements[globalDofIdx] = remainderDisplacements[remainderDofIdx];
                }

                // For boundary dofs we take average across subdomains with respect to multiplicity or stiffness. 
                double[] boundaryDofCoeffs = distribution.CalcBoundaryDofCoefficients(subdomain);
                for (int i = 0; i < dofSeparator.BoundaryDofIndices[id].Length; ++i)
                {
                    int remainderDofIdx = dofSeparator.BoundaryDofIndices[id][i];
                    int globalDofIdx = freeToGlobalDofs[remainderToFreeDofs[remainderDofIdx]];
                    globalDisplacements[globalDofIdx] += remainderDisplacements[remainderDofIdx] * boundaryDofCoeffs[i];
                }
            }

            // Corner dofs are copied without averaging.
            for (int i = 0; i < dofSeparator.NumGlobalCornerDofs; ++i)
            {
                int globalDofIdx = dofSeparator.GlobalCornerToFreeDofMap[i];
                globalDisplacements[globalDofIdx] = globalCornerDisplacements[i];
            }

            return globalDisplacements;
        }

        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainDisplacements)
        {
            var globalDisplacements = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);

            // Remainder dofs
            foreach (var subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                int[] remainderToSubdomainDofs = dofSeparator.RemainderDofIndices[id];
                int[] cornerToSubdomainDofs = dofSeparator.CornerDofIndices[id];
                IVectorView freeDisplacements = subdomainDisplacements[id]; //TODO: benchmark the performance if this was concrete Vector

                // Internal dofs: We copy them without averaging.
                foreach (int remainderDofIdx in dofSeparator.InternalDofIndices[id])
                {
                    int subdomainDofIdx = remainderToSubdomainDofs[remainderDofIdx];
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                    globalDisplacements[globalDofIdx] = freeDisplacements[subdomainDofIdx];
                }

                // Boundary remainder dofs: We take average across subdomains with respect to multiplicity or stiffness. 
                double[] boundaryDofCoeffs = distribution.CalcBoundaryDofCoefficients(subdomain);
                for (int i = 0; i < dofSeparator.BoundaryDofIndices[id].Length; ++i)
                {
                    int subdomainDofIdx = remainderToSubdomainDofs[dofSeparator.BoundaryDofIndices[id][i]];
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                    globalDisplacements[globalDofIdx] += freeDisplacements[subdomainDofIdx] * boundaryDofCoeffs[i];
                }

                // Boundary corner dofs: We copy without averaging.
                foreach (int subdomainDofIdx in cornerToSubdomainDofs)
                {
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];

                    // This will overwrite the value from a previous subdomain, but these values are the same.
                    globalDisplacements[globalDofIdx] = freeDisplacements[subdomainDofIdx]; 
                }
            }

            return globalDisplacements;
        }

        public Vector GatherGlobalForces(Dictionary<int, IVectorView> subdomainForces)
        {
            var globalForces = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                int[] subdomainFreeToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                IVectorView forces = subdomainForces[id]; //TODO: benchmark the performance if this was concrete Vector

                for (int i = 0; i < forces.Length; ++i)
                {
                    // Internal forces will be copied (which is identical to adding 0 + single value).
                    // Boundary remainder forces will be summed. Previously we had distributed them depending on 
                    // homogeneity / heterogeneity (e.g. Ftot = 0.4 * Ftot + 0.6 * Ftot) and now we sum them. 
                    // Boundary corner forces are also summed. Previously we had also distributed them equally irregardless of 
                    // homogeneity / heterogeneity (e.g. Ftot = 0.5 * Ftot + 0.5 * Ftot) and now we sum them.
                    globalForces[subdomainFreeToGlobalDofs[i]] += forces[i];
                }
            }
            return globalForces;
        }
    }
}
