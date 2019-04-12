using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IReferenceVolumeElement
    {
        void ApplyBoundaryConditions();
        IMatrixView CalculateKinematicRelationsMatrix(ISubdomain subdomain);
        double CalculateRveVolume();
    }
}
