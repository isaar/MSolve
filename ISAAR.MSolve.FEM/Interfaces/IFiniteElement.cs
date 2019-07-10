using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IFiniteElement : IElementType
    {
        int ID { get; }
        ElementDimensions ElementDimensions { get; }
        bool MaterialModified { get; }
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
