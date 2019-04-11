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

		public AsymmetricDofOrderer(IAsymmetricDofOrderingStrategy rowOrderingStrategy)
		{
			_rowOrderingStrategy = rowOrderingStrategy;
		}

		public IGlobalFreeDofOrdering OrderDofs(IStructuralAsymmetricModel model)
		{
			var subdomain = model.Subdomains.First();
			(int numSubdomainFreeRowDofs, DofTable subdomainFreeRowDofs) = _rowOrderingStrategy.OrderSubdomainDofs(subdomain);
			ISubdomainFreeDofOrdering subdomainRowOrdering= new SubdomainFreeRowDofOrderingGeneral(numSubdomainFreeRowDofs, subdomainFreeRowDofs);

			return  new GlobalFreeDofOrderingSingle((ISubdomain_v2)subdomain, subdomainRowOrdering);
		}
	}
}
