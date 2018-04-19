using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public enum ElementDimensions
    {
        Unknown = 0,
        OneD = 1,
        TwoD = 2,
        ThreeD = 3
    }

    public interface IFiniteElement //TODOMaria elemental loads should be calculated like this: the one who calls CalculateStresses et.c. should provide nodal displacements compliant to the constraints
    {
        int ID { get; }
        ElementDimensions ElementDimensions { get; }
        IFiniteElementDOFEnumerator DOFEnumerator { get; set; }
        IList<IList<DOFType>> GetElementDOFTypes(Element element);
        bool MaterialModified { get; }
        IMatrix2D StiffnessMatrix(Element element);
        IMatrix2D MassMatrix(Element element);
        IMatrix2D DampingMatrix(Element element);
        void ResetMaterialModified();
        Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForcesForLogging(Element element, double[] localDisplacements);
        double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads);
        void SaveMaterialState();
        void ClearMaterialState();

        void ClearMaterialStresses();
    }
}
