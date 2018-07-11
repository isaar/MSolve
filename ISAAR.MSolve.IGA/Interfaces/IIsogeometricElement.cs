using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.IGA.Interfaces
{

    public interface IIsogeometricElement: IElementType
	{
        int ID { get; }
        ElementDimensions ElementDimensions { get; }
	    IElementDOFEnumerator DOFEnumerator { get; set; }
        //IList<IList<DOFType>> GetElementDOFTypes(IElement element);
        bool MaterialModified { get; }
        //IMatrix2D StiffnessMatrix(Element element);
        //IMatrix2D MassMatrix(Element element);
        //IMatrix2D DampingMatrix(Element element);
        Dictionary<int, double> CalculateLoadingCondition(Element element,Edge edge, NeumannBoundaryCondition neumann);
        Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann);
        Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure);
        Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure);
        void ResetMaterialModified();
        Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForcesForLogging(Element element, double[] localDisplacements);
        //double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads);
        //void SaveMaterialState();
        void ClearMaterialState();

        //void ClearMaterialStresses();
    }
}
