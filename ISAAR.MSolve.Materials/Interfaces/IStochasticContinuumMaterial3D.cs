using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IStochasticContinuumMaterial3D : IContinuumMaterial3D
    {
        IStochasticMaterialCoefficientsProvider CoefficientsProvider { get; set; }
        IMatrix2D GetConstitutiveMatrix(double[] coordinates);
    }
}
