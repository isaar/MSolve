using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
	public class GlobalFreeDofAsymmetricOrderingSingle : IGlobalFreeDofOrdering
	{
		private readonly ISubdomain _subdomain;
		private readonly int[] _subdomainToGlobalDofMap;

		public GlobalFreeDofAsymmetricOrderingSingle(ISubdomain subdomain,
			ISubdomainFreeDofOrdering subdomainRowOrdering, ISubdomainFreeDofOrdering subdomainColOrdering)
		{
			_subdomain = subdomain;
			this.NumGlobalFreeDofs = subdomainColOrdering.NumFreeDofs;
			this.GlobalFreeDofs = subdomainColOrdering.FreeDofs;

		}


		public DofTable GlobalFreeDofs { get; }
		public int NumGlobalFreeDofs { get; }
		public IReadOnlyDictionary<ISubdomain, ISubdomainFreeDofOrdering> SubdomainDofOrderings { get; }

		public void AddVectorSubdomainToGlobal(ISubdomain subdomain, IVectorView subdomainVector,
			IVector globalVector)
		{
			throw new NotImplementedException();
		}

		public void AddVectorSubdomainToGlobalMeanValue(ISubdomain subdomain, IVectorView subdomainVector,
			IVector globalVector)
		{
			throw new NotImplementedException();
		}

		public void ExtractVectorSubdomainFromGlobal(ISubdomain subdomain, IVectorView globalVector,
			IVector subdomainVector)
		{
			throw new NotImplementedException();
		}

		public int[] MapFreeDofsSubdomainToGlobal(ISubdomain subdomain)
		{
			throw new NotImplementedException();
		}
	}
}