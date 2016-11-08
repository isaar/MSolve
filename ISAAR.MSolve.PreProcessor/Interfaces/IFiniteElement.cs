using System;
using System.Collections.Generic;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public enum ElementDimensions
    {
        Unknown = 0,
        OneD = 1,
        TwoD = 2,
        ThreeD = 3
    }

    public interface IFiniteElement
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
