using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Analyzers.Multiscale
{
    public interface IReferenceVolumeElement
    {
        void ApplyBoundaryConditions(IStructuralModel_v2 model);
        double CalculateRveVolume();
    }
}
