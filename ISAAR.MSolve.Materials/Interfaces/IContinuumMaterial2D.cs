using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IContinuumMaterial2D : IFiniteElementMaterial
    {
        double[] Stresses { get; }
        IMatrixView ConstitutiveMatrix { get; }
        void UpdateMaterial(double[] strains);
    }
}
