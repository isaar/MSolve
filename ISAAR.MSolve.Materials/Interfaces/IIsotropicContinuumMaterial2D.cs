using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IIsotropicContinuumMaterial2D : IContinuumMaterial2D
    {
        double YoungModulus { get; set; }
        double PoissonRatio { get; set; }
    }
}
