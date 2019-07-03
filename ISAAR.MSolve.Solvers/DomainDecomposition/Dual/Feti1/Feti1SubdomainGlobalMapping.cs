using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

//TODO: the operations described here are expressed as matrix-vector multiplications with the projection & scaling matrices 
//      called Lpb in FETI bibliography. Should I also express them similarly? In that case dedicated classes are needed for 
//      the Lpb matrices, but the code for the conversions is the same for homogeneous/heterogeneous. Only the creation of the
//      Lpb matrices depends on the stiffness distribution, which is expected and needed elsewhere too.
//TODO: These methods overlap with IGlobalDofOrdering. Perhaps they should be moved there, while the implementing classes remain 
//      here as strategies.
//TODO: Should the vectors be Vector instead of IVectorView?
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1
{
    public class Feti1SubdomainGlobalMapping
    {
        private readonly IStiffnessDistribution distribution;
        private readonly Feti1DofSeparator dofSeparator;
        private readonly IStructuralModel model;

        public Feti1SubdomainGlobalMapping(IStructuralModel model, Feti1DofSeparator dofSeparator,
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

            var vectors = new Dictionary<int, SparseVector>();
            foreach (var idSubdomains in subdomains)
            {
                int id = idSubdomains.Key;
                int numSubdomainDofs = idSubdomains.Value.FreeDofOrdering.NumFreeDofs;
                vectors[id] = SparseVector.CreateFromDictionary(numSubdomainDofs, subdomainLoads[id]);
            }
            return vectors;
        }

        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainDisplacements)
        {
            var globalDisplacements = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            foreach (var subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                IVectorView displacements = subdomainDisplacements[id]; //TODO: benchmark the performance if this was concrete Vector

                // Internal dofs are copied without averaging.
                foreach (int internalDof in dofSeparator.InternalDofIndices[id])
                {
                    int globalDofIdx = subdomainToGlobalDofs[internalDof];
                    globalDisplacements[globalDofIdx] = displacements[internalDof];
                }

                // For boundary dofs we take average across subdomains with respect to multiplicity or stiffness. 
                double[] boundaryDofCoeffs = distribution.CalcBoundaryDofCoefficients(subdomain);
                for (int i = 0; i < dofSeparator.BoundaryDofIndices[id].Length; ++i)
                {
                    int subdomainDofIdx = dofSeparator.BoundaryDofIndices[id][i];
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                    globalDisplacements[globalDofIdx] += displacements[subdomainDofIdx] * boundaryDofCoeffs[i];
                }
            }
            return globalDisplacements;
        }

        public Vector GatherGlobalForces(Dictionary<int, IVectorView> subdomainForces)
        {
            // For boundary dofs we sum the contributions from each subdomain. There is no averaging here, unlike the
            // conversion of nodal loads from global to subdomain dofs. That averaging only holds for nodal loads, while
            // loads calculated from elements (e.g. Dirichlet, Neumann) are applied directly to the subdomain to which the 
            // element belongs. To find the global forces, we just need to sum the contributions from each subdomain. This
            // is also correct for the conversion of nodal loads from subdomain to global dofs, since they have been already
            // distributed to each subdomain, according to the multiplicity or relative stiffness of the corresponding dof.

            var globalForces = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            foreach (var subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                IVectorView forces = subdomainForces[id]; //TODO: benchmark the performance if this was concrete Vector

                for (int i = 0; i < forces.Length; ++i)
                {
                    // Boundary forces will be summed. 
                    // Internal forces will be copied (which is identical to adding 0 + single value).
                    globalForces[subdomainToGlobalDofs[i]] += forces[i];
                }
            }
            return globalForces;
        }

        #region incorrect implementation
        // This is version is more efficient, because it only transfers the sum (double) per subdomain, instead of a vector 
        // corresponding to the boundary dofs of the whole model. However it is correct if there are only nodal loads, since it 
        // assumes that the force at each boundary dof originated as a nodal loads, that was distributed (according to the 
        // multiplicity for homogeneous problems or the relative stiffness for heterogeneous). This doesn't hold for forces
        // resulting from loads applied to elements (e.g. Dirichlet, Neumann), because they need to be summed without any 
        // averaging.
        //TODO: Could this logic be used as an optimization for the case it is correct?
        //private double CalculateGlobalForcesNorm_INCORRECT(Dictionary<int, IVectorView> subdomainForces)
        //{
        //    //TODO: this should be used for non linear analyzers as well (instead of building the global RHS)
        //    //TODO: Is this correct? For the residual, it would be wrong to find f-K*u for each subdomain and then call this.

        //    double globalSum = 0.0;
        //    foreach (ISubdomain subdomain in model.Subdomains)
        //    {
        //        int id = subdomain.ID;
        //        double subdomainSum = 0.0;
        //        IVectorView forces = subdomainForces[id];
        //        foreach (int internalDof in dofSeparator.InternalDofIndices[id]) subdomainSum += forces[internalDof];
        //        for (int i = 0; i < dofSeparator.BoundaryDofIndices[id].Length; ++i)
        //        {
        //            // E.g. 2 subdomains. Originally: sum += f * f. Now: sum += (2*f/2)(2*f/2)/2 + (2*f/2)(2*f/2)/2
        //            // WARNING: This works only if nodal loads are distributed evenly among subdomains.
        //            int multiplicity = dofSeparator.BoundaryDofMultiplicities[id][i];
        //            double totalForce = forces[dofSeparator.BoundaryDofIndices[id][i]] * multiplicity;
        //            subdomainSum += (totalForce * totalForce) / multiplicity;
        //        }
        //        globalSum += subdomainSum;
        //    }
        //    return Math.Sqrt(globalSum);
        //}
        #endregion
    }
}
