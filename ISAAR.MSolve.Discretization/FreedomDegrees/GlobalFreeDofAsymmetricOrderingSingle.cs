using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
	public class GlobalFreeDofAsymmetricOrderingSingle : IGlobalFreeDofOrdering
	{
		private readonly ISubdomain_v2 _subdomain;
		private readonly int[] _subdomainToGlobalDofMap;

		public GlobalFreeDofAsymmetricOrderingSingle(ISubdomain_v2 subdomain,
			ISubdomainFreeDofOrdering subdomainRowOrdering, ISubdomainFreeDofOrdering subdomainColOrdering)
		{
			_subdomain = subdomain;
			this.NumGlobalFreeDofs = subdomainColOrdering.NumFreeDofs;
			this.GlobalFreeDofs = subdomainColOrdering.FreeDofs;

		}


		public DofTable GlobalFreeDofs { get; }
		public int NumGlobalFreeDofs { get; }
		public IReadOnlyDictionary<ISubdomain_v2, ISubdomainFreeDofOrdering> SubdomainDofOrderings { get; }

		public void AddVectorSubdomainToGlobal(ISubdomain_v2 subdomain, IVectorView subdomainVector,
			IVector globalVector)
		{
			throw new NotImplementedException();
		}

		public void AddVectorSubdomainToGlobalMeanValue(ISubdomain_v2 subdomain, IVectorView subdomainVector,
			IVector globalVector)
		{
			throw new NotImplementedException();
		}

		public void ExtractVectorSubdomainFromGlobal(ISubdomain_v2 subdomain, IVectorView globalVector,
			IVector subdomainVector)
		{
			throw new NotImplementedException();
		}

		public int[] MapFreeDofsSubdomainToGlobal(ISubdomain_v2 subdomain)
		{
			throw new NotImplementedException();
		}
	}
}