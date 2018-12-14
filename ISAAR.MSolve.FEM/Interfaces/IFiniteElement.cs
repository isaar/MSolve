using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IFiniteElement : IElementType
    {
        int ID { get; }
        ElementDimensions ElementDimensions { get; }
        //IElementDOFEnumerator DOFEnumerator { get; set; }
        //IList<IList<DOFType>> GetElementDOFTypes(Element element);
        bool MaterialModified { get; }
        //IMatrix2D StiffnessMatrix(Element element);
        //IMatrix2D MassMatrix(Element element);
        //IMatrix2D DampingMatrix(Element element);
        void ResetMaterialModified();
        Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForcesForLogging(Element element, double[] localDisplacements);
        double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads);
        void SaveMaterialState();
        void ClearMaterialState();

        void ClearMaterialStresses(); //TODO this is only for structural problems.
    }
}
