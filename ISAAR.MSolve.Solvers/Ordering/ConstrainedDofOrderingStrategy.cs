using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: This could be done simultaneously with ordering the free dofs, to improve performance.
namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Orders the constrained dofs of a subdomain, independendtly from the free ones.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class ConstrainedDofOrderingStrategy
    {
        /// <summary>
        /// Orders the constrained freedom degrees of one of the model's subdomains.
        /// </summary>
        /// <param name="subdomain">A subdomain of the whole model.</param>
        internal (int numSubdomainConstrainedDofs, DofTable subdomainConstrainedDofs) OrderSubdomainDofs(ISubdomain_v2 subdomain)
        {
            var constrainedDofs = new DofTable();
            int dofCounter = 0;
            foreach (INode node in subdomain.Nodes)
            {
                foreach (Constraint constraint in node.Constraints) constrainedDofs[node, constraint.DOF] = dofCounter++;
            }
            return (dofCounter, constrainedDofs);
        }
    }
}
