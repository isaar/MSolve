using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1
{
    /// <summary>
    /// Calculates the signed boolean matrices of the equations that enforce continuity between the multiple instances of 
    /// boundary dofs.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Feti1LagrangeMultipliersEnumerator : LagrangeMultipliersEnumeratorBase
    {
        internal Feti1LagrangeMultipliersEnumerator(ICrosspointStrategy crosspointStrategy, Feti1DofSeparator dofSeparator) 
            : base(crosspointStrategy, dofSeparator)
        {
        }

        /// <summary>
        /// Only <see cref="BooleanMatrices"/> will be explicitly created. <see cref="LagrangeMultipliers"/> will not.
        /// For use in homogeneous problems, where we do not need that much info about lagrange multipliers and boundary dofs.
        /// </summary>
        /// <param name="model"></param>
        public void DefineBooleanMatrices(IStructuralModel model)
        {
            var numFreeDofs = new Dictionary<int, int>();
            var freeDofOrderings = new Dictionary<int, DofTable>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                numFreeDofs[subdomain.ID] = subdomain.FreeDofOrdering.NumFreeDofs;
                freeDofOrderings[subdomain.ID] = subdomain.FreeDofOrdering.FreeDofs;
            }
            base.DefineBooleanMatrices(model, numFreeDofs, freeDofOrderings);
        }

        /// <summary>
        /// Creates both <see cref="BooleanMatrices"/> and <see cref="LagrangeMultipliers"/>. For use in heterogeneous problems.
        /// </summary>
        /// <param name="model"></param>
        public void DefineLagrangesAndBooleanMatrices(IStructuralModel model)
        {
            var numFreeDofs = new Dictionary<int, int>();
            var freeDofOrderings = new Dictionary<int, DofTable>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                numFreeDofs[subdomain.ID] = subdomain.FreeDofOrdering.NumFreeDofs;
                freeDofOrderings[subdomain.ID] = subdomain.FreeDofOrdering.FreeDofs;
            }
            base.DefineLagrangesAndBooleanMatrices(model, numFreeDofs, freeDofOrderings);
        }
    }
}
