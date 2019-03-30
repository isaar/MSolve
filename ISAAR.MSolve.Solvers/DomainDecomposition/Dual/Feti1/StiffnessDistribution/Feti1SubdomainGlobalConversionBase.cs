using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: Also abstract DistributeNodalLoads as much as possible. Resuing the LumpedBoundaryStiffnesses of heterogeneous 
//      distribution might help. Ideally this class would not be abstract. Instead it should be injected with the relative 
//      stiffness / multiplicity for each dof from the IStiffnessDistribution and just use that.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution
{
    internal abstract class Feti1SubdomainGlobalConversionBase : ISubdomainGlobalConversion
    {
        protected readonly IStructuralModel_v2 model;
        protected readonly Feti1DofSeparator dofSeparator;

        protected Feti1SubdomainGlobalConversionBase(IStructuralModel_v2 model, Feti1DofSeparator dofSeparator)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
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

        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainDisplacements)
        {
            var globalDisplacements = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            foreach (var subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                IVectorView displacements = subdomainDisplacements[id]; //TODO: benchmark the performance if this was concrete Vector

                // Internal dofs are copied as is.
                foreach (int internalDof in dofSeparator.InternalDofIndices[id])
                {
                    int globalDofIdx = subdomainToGlobalDofs[internalDof];
                    globalDisplacements[globalDofIdx] = displacements[internalDof];
                }

                // For boundary dofs we take average across subdomains with respect to multiplicity or stiffness. 
                double[] boundaryDofCoeffs = CalcBoundaryDofMultipliers(subdomain);
                for (int i = 0; i < dofSeparator.BoundaryDofIndices[id].Length; ++i)
                {
                    int multiplicity = dofSeparator.BoundaryDofMultiplicities[id][i];
                    int subdomainDofIdx = dofSeparator.BoundaryDofIndices[id][i];
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                    globalDisplacements[globalDofIdx] += displacements[subdomainDofIdx] * boundaryDofCoeffs[i];
                }
            }
            return globalDisplacements;
        }

        public Vector GatherGlobalForces(Dictionary<int, IVectorView> subdomainForces)
        {
            var globalForces = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            foreach (var subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                IVectorView forces = subdomainForces[id]; //TODO: benchmark the performance if this was concrete Vector

                // Internal dofs are copied as is.
                foreach (int internalDof in dofSeparator.InternalDofIndices[id])
                {
                    int globalDofIdx = subdomainToGlobalDofs[internalDof];
                    globalForces[globalDofIdx] = forces[internalDof];
                }

                // For boundary dofs we sum the contributions from each subdomain. There is no averaging here, unlike the
                // conversion of nodal loads from global to subdomain dofs. That averaging only holds for nodal loads, while
                // loads calculated from elements (e.g. Dirichlet, Neumann) are applied directly to the subdomain to which the 
                // element belongs. To find the global forces, we just need to sum the contributions from each subdomain. This
                // is also correct for the conversion of nodal loads from subdomain to global dofs, since they have been already
                // distributed to each subdomain, according to the multiplicity or relative stiffness of the corresponding dof.
                for (int i = 0; i < dofSeparator.BoundaryDofIndices[id].Length; ++i)
                {
                    int subdomainDofIdx = dofSeparator.BoundaryDofIndices[id][i];
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                    globalForces[globalDofIdx] += forces[subdomainDofIdx];
                }
            }
            return globalForces;
        }

        protected Dictionary<int, SparseVector> BuildForceVectors(IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems,
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

        public abstract Dictionary<int, SparseVector> DistributeNodalLoads(
            IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems, Table<INode, DOFType, double> globalNodalLoads);
        protected abstract double[] CalcBoundaryDofMultipliers(ISubdomain_v2 subdomain);

        #region incorrect implementation
        // This is version is more efficient, because it only transfers the sum (double) per subdomain, instead of a vector 
        // corresponding to the boundary dofs of the whole model. However it is correct if there are only nodal loads, since it 
        // assumes that the force at each boundary dof originated as a nodal loads, that was distributed (according to the 
        // multiplicity for homogeneous problems or the relative stiffness for heterogeneous). This doesn't hold for forces
        // resulting from loads applied to elements (e.g. Dirichlet, Neumann), because they need to be summed without any 
        // averaging.
        //TODO: Could this logic be used as an optimization for the case it is correct?
        private double CalculateGlobalForcesNorm_INCORRECT(Dictionary<int, IVectorView> subdomainForces)
        {
            //TODO: this should be used for non linear analyzers as well (instead of building the global RHS)
            //TODO: Is this correct? For the residual, it would be wrong to find f-K*u for each subdomain and then call this.

            double globalSum = 0.0;
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                double subdomainSum = 0.0;
                IVectorView forces = subdomainForces[id];
                foreach (int internalDof in dofSeparator.InternalDofIndices[id]) subdomainSum += forces[internalDof];
                for (int i = 0; i < dofSeparator.BoundaryDofIndices[id].Length; ++i)
                {
                    // E.g. 2 subdomains. Originally: sum += f * f. Now: sum += (2*f/2)(2*f/2)/2 + (2*f/2)(2*f/2)/2
                    // WARNING: This works only if nodal loads are distributed evenly among subdomains.
                    int multiplicity = dofSeparator.BoundaryDofMultiplicities[id][i];
                    double totalForce = forces[dofSeparator.BoundaryDofIndices[id][i]] * multiplicity;
                    subdomainSum += (totalForce * totalForce) / multiplicity;
                }
                globalSum += subdomainSum;
            }
            return Math.Sqrt(globalSum);
        }
        #endregion
    }
}
