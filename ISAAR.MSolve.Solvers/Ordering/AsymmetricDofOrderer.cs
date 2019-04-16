using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Assemblers.Collocation;

namespace ISAAR.MSolve.Solvers.Ordering
{
	/// <summary>
	/// Orders the unconstrained freedom degrees of each subdomain and the whole model.
	/// It is used for ordering elements that produce non-symmetric matrices or with different row-column ordering.
	/// </summary>
	
	public class AsymmetricDofOrderer: IAsymmetricDofOrderer
    {
		private readonly IAsymmetricDofOrderingStrategy _rowOrderingStrategy;
        private readonly ConstrainedDofOrderingStrategy _constrainedOrderingStrategy;
        private readonly bool _cacheElementToSubdomainDofMaps = true;

        public AsymmetricDofOrderer(IAsymmetricDofOrderingStrategy rowOrderingStrategy, bool cacheElementToSubdomainMaps=true)
        {
            _cacheElementToSubdomainDofMaps = cacheElementToSubdomainMaps;
            _constrainedOrderingStrategy= new ConstrainedDofOrderingStrategy();
			_rowOrderingStrategy = rowOrderingStrategy;
		}

		public IGlobalFreeDofOrdering OrderFreeDofs(IStructuralAsymmetricModel model)
		{
			var subdomain = model.Subdomains.First();
			(int numSubdomainFreeRowDofs, DofTable subdomainFreeRowDofs) = _rowOrderingStrategy.OrderSubdomainDofs(subdomain);
			ISubdomainFreeDofOrdering subdomainRowOrdering= new SubdomainFreeRowDofOrderingGeneral(numSubdomainFreeRowDofs, subdomainFreeRowDofs);

			return  new GlobalFreeDofOrderingSingle((ISubdomain)subdomain, subdomainRowOrdering);
		}

        public ISubdomainConstrainedDofOrdering OrderConstrainedDofs(ISubdomain subdomain)
        {
            (int numConstrainedDofs, DofTable constrainedDofs) =
                _constrainedOrderingStrategy.OrderSubdomainDofs(subdomain);
            if (_cacheElementToSubdomainDofMaps)
            {
                return new SubdomainConstrainedDofOrderingCaching(numConstrainedDofs, constrainedDofs);
            }
            return new SubdomainConstrainedDofOrderingGeneral(numConstrainedDofs,constrainedDofs);
        }
	}
}
