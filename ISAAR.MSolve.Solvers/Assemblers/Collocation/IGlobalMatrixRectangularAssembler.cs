using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
	public interface IGlobalMatrixRectangularAssembler<TMatrix> where TMatrix:IMatrix
	{
		TMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofRowOrdering, ISubdomainFreeDofOrdering dofColOrdering,IEnumerable<IElement> elements,
			IElementMatrixProvider matrixProvider);

        (TMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree, IMatrixView
            matrixConstrConstr)
            BuildGlobalSubmatrices(
                ISubdomainFreeDofOrdering freeDofRowOrdering,
                ISubdomainFreeDofOrdering freeDofColOrdering,
                ISubdomainConstrainedDofOrdering constrainedDofRowOrdering,
                ISubdomainConstrainedDofOrdering constrainedDofColOrdering,
                IEnumerable<IElement> elements, IElementMatrixProvider matrixProvider);

        void HandleDofOrderingWillBeModified();
	}
}