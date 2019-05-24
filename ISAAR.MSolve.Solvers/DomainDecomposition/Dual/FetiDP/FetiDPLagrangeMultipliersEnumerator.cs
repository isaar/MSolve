using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: This works for other dual solvers, not only FETI-1
//TODO: This should only calculate them. Another object should manage them.
//TODO: The enumation code is quite complex and error prone. It should be simplified and decomposed into smaller methods, as
//      much as possible without sacrificing performance.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    /// <summary>
    /// Calculates the signed boolean matrices of the equations that enforce continuity between the multiple instances of 
    /// boundary dofs.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class FetiDPLagrangeMultipliersEnumerator : LagrangeMultipliersEnumeratorBase
    {
        private readonly FetiDPDofSeparator dofSeparator;

        internal FetiDPLagrangeMultipliersEnumerator(ICrosspointStrategy crosspointStrategy, FetiDPDofSeparator dofSeparator) 
            : base(crosspointStrategy, dofSeparator)
        {
            this.dofSeparator = dofSeparator;
        }

        /// <summary>
        /// Only <see cref="BooleanMatrices"/> will be explicitly created. <see cref="LagrangeMultipliers"/> will not.
        /// For use in homogeneous problems, where we do not need that much info about lagrange multipliers and boundary dofs.
        /// </summary>
        /// <param name="model"></param>
        public void DefineBooleanMatrices(IStructuralModel model)
        {
            var numRemainderDofs = new Dictionary<int, int>();
            var remainderDofOrderings = new Dictionary<int, DofTable>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                numRemainderDofs[subdomain.ID] = dofSeparator.RemainderDofIndices[subdomain.ID].Length;
                remainderDofOrderings[subdomain.ID] = dofSeparator.RemainderDofOrderings[subdomain.ID];
            }
            base.DefineBooleanMatrices(model, numRemainderDofs, remainderDofOrderings);
        }

        /// <summary>
        /// Creates both <see cref="BooleanMatrices"/> and <see cref="LagrangeMultipliers"/>. For use in heterogeneous problems.
        /// </summary>
        /// <param name="model"></param>
        public void DefineLagrangesAndBooleanMatrices(IStructuralModel model)
        {
            var numRemainderDofs = new Dictionary<int, int>();
            var remainderDofOrderings = new Dictionary<int, DofTable>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                numRemainderDofs[subdomain.ID] = dofSeparator.RemainderDofIndices[subdomain.ID].Length;
                remainderDofOrderings[subdomain.ID] = dofSeparator.RemainderDofOrderings[subdomain.ID];
            }
            base.DefineLagrangesAndBooleanMatrices(model, numRemainderDofs, remainderDofOrderings);
        }
    }
}
