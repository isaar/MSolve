using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: The semantics of boundary & internal dofs is different for each domain decomposition method. Perhaps this interface 
//      should not exist. Instead the Feti1DofSeparator, FetiDPDofSeparator, etc. would follow the same design and compose the 
//      same utility classes. Then the classes that actually know the semantics of boundary & internal dofs of each FETI method
//      (e.g. Feti1Solver, Feti1PreconditionerFactory) will have access to concrete DofSeparators, access appropriate boundary 
//      dof data from there and pass it to classes that work with that data, instead of IDofSeparator 
//      (e.g. Homogeneous/HeterogeneousStiffnessDistribution).
namespace ISAAR.MSolve.Solvers.DomainDecomposition
{
    public interface IDofSeparator
    {
        /// <summary>
        /// The indices of dofs, where Lagrange multipliers will be applied, into a sequence of subdomain dofs. Depending on the 
        /// domain decomposition method, that sequence could be all free boundary dofs (e.g. FETI-1) or a subset of them 
        /// (e.g. boundary dofs minus corner dofs in FETI-DP).
        /// </summary>
        Dictionary<int, int[]> BoundaryDofIndices { get; }

        /// <summary>
        /// Dofs of each subdomain where Lagrange multipliers will be applied. Depending on the domain decomposition method, 
        /// these could be all free boundary dofs (e.g. FETI-1) or a subset of them (e.g. boundary dofs minus corner dofs in 
        /// FETI-DP).
        /// </summary>
        Dictionary<int, (INode node, IDofType dofType)[]> BoundaryDofs { get; }

        /// <summary>
        /// Dofs where Lagrange multipliers will be applied. Depending on the domain decomposition method, these could be all
        /// free boundary dofs (e.g. FETI-1) or a subset of them (e.g. boundary dofs minus corner dofs in FETI-DP).
        /// </summary>
        Dictionary<INode, IDofType[]> GlobalBoundaryDofs { get; }

        /// <summary>
        /// The indices of dofs, which only occur for each subdomain, into a sequence of subdomain dofs. Depending on the 
        /// domain decomposition method, that sequence could be all free boundary dofs (e.g. FETI-1) or a subset of them 
        /// (e.g. boundary dofs minus corner dofs in FETI-DP).
        /// </summary>
        Dictionary<int, int[]> InternalDofIndices { get; }
    }
}
