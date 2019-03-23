using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IFiniteElement_v2 : IElementType_v2
    {
        int ID { get; }
        ElementDimensions ElementDimensions { get; }
        bool MaterialModified { get; }
        void ResetMaterialModified();
        Tuple<double[], double[]> CalculateStresses(Element_v2 element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForces(Element_v2 element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForcesForLogging(Element_v2 element, double[] localDisplacements);
        double[] CalculateAccelerationForces(Element_v2 element, IList<MassAccelerationLoad> loads);
        void SaveMaterialState();
        void ClearMaterialState();

        void ClearMaterialStresses(); //TODO this is only for structural problems.
    }
}
