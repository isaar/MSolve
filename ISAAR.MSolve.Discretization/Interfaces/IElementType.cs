using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
	public enum ElementDimensions
	{
		Unknown = 0,
		OneD = 1,
		TwoD = 2,
		ThreeD = 3
	}

	public interface IElementType
    {
		IElementDOFEnumerator DOFEnumerator { get; set; }
	    IMatrix2D StiffnessMatrix(IElement element);
	    IMatrix2D MassMatrix(IElement element);
	    IMatrix2D DampingMatrix(IElement element);
	    IList<IList<DOFType>> GetElementDOFTypes(IElement element);
	}
}
