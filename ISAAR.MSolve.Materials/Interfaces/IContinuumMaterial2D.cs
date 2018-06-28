using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IContinuumMaterial2D : IFiniteElementMaterial
    {
        StressStrainVectorContinuum2D Stresses { get; }
        ElasticityTensorContinuum2D ConstitutiveMatrix { get; }
        void UpdateMaterial(StressStrainVectorContinuum2D strains);
        void ClearState();
        void SaveState();
        void ClearStresses();
    }
}
