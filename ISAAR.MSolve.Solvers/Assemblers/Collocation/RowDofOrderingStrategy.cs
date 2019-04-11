using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
	public class RowDofOrderingStrategy : IAsymmetricDofOrderingStrategy
	{
		private readonly IReadOnlyList<DOFType> dofsPerNode= new List<DOFType> {DOFType.X, DOFType.Y };

		public (int numGlobalFreeDofs, DofTable globalFreeDofs) OrderGlobalDofs(IStructuralAsymmetricModel model)
			=> OrderFreeDofsOfElementSet(model.Elements, model.Constraints);

        private (int numGlobalFreeDofs, DofTable globalFreeDofs) OrderFreeDofsOfElementSet(IReadOnlyList<IElement_v2> elements, Table<INode, DOFType, double> constraints)
        {
            var freeDofs = new DofTable();
            int dofCounter = 0;
            foreach (var element in elements)
            {
                var collocationPoint = ((ICollocationElement)element).CollocationPoint;
                bool isNodeConstrained = collocationPoint.Constraints?.Count != 0;
                foreach (DOFType dof in dofsPerNode)
                {
                    bool isDofConstrained = isNodeConstrained && collocationPoint.Constraints!=null&& collocationPoint.Constraints.Any(c=>c.DOF==dof);
                    if (!isDofConstrained)
                        freeDofs[collocationPoint, dof] = dofCounter++;
                }
            }
            return (dofCounter, freeDofs);
        }


		public (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) OrderSubdomainDofs(IAsymmetricSubdomain subdomain)
			=> OrderFreeDofsOfElementSet(((ISubdomain_v2)subdomain).Elements, ((ISubdomain_v2)subdomain).Constraints);
	}
}