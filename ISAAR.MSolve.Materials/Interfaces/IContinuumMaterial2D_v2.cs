using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IContinuumMaterial2D_v2 : IFiniteElementMaterial
    {
        IVectorView Stresses { get; }
        IMatrixView ConstitutiveMatrix { get; }
        void UpdateMaterial(IVectorView strains);
    }
}
