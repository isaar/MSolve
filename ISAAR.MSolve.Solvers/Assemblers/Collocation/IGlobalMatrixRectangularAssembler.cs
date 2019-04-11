using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
	public interface IGlobalMatrixRectangularAssembler<TMatrix> where TMatrix:IMatrix
	{
		TMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofRowOrdering, ISubdomainFreeDofOrdering dofColOrdering,IEnumerable<IElement_v2> elements,
			IElementMatrixProvider_v2 matrixProvider);

		void HandleDofOrderingWillBeModified();
	}
}