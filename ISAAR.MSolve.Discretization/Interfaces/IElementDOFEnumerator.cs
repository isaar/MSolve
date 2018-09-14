using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
	public interface IElementDOFEnumerator
	{
		IList<IList<DOFType>> GetDOFTypes(IElement element);
		IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(IElement element);
		IList<INode> GetNodesForMatrixAssembly(IElement element);
		IMatrix2D GetTransformedMatrix(IMatrix2D matrix);
		double[] GetTransformedVector(double[] vector);
	}
}