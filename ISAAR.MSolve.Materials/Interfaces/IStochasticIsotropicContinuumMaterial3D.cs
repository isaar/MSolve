using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IStochasticIsotropicContinuumMaterial3D : IIsotropicContinuumMaterial3D //maybe we don't need it?? the interface members are unused
    {
        IStochasticMaterialCoefficientsProvider CoefficientsProvider { get; set; }
        ElasticityTensorContinuum3D GetConstitutiveMatrix(double[] coordinates);
    }
}
