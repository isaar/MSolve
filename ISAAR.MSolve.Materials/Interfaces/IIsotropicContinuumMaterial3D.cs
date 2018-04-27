using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IIsotropicContinuumMaterial3D : IContinuumMaterial3D
    {
        double YoungModulus { get; set; }
        double PoissonRatio { get; set; }
    }
}
