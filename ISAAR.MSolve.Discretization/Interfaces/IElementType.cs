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
	    int ID { get; }
	    ElementDimensions ElementDimensions { get; }
		IElementDOFEnumerator DOFEnumerator { get; set; }
		IList<IList<DOFType>> GetElementDOFTypes(IElement element);
	    bool MaterialModified { get; }
	    IMatrix2D StiffnessMatrix(IElement element);
	    IMatrix2D MassMatrix(IElement element);
	    IMatrix2D DampingMatrix(IElement element);
	    void ResetMaterialModified();
	    Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements);
	    double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements);
	    double[] CalculateForcesForLogging(IElement element, double[] localDisplacements);
	    void SaveMaterialState();
	    void ClearMaterialState();
	    void ClearMaterialStresses();
	}
}
