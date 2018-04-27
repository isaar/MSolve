
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IContinuumMaterial3D : IFiniteElementMaterial
    {
        StressStrainVectorContinuum3D Stresses { get; }
        ElasticityTensorContinuum3D ConstitutiveMatrix { get; }
        void UpdateMaterial(StressStrainVectorContinuum3D strains);
        void ClearState();
        void SaveState();
        void ClearStresses();
    }
}
