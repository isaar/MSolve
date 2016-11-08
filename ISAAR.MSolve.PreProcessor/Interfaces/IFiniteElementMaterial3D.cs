using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public interface IFiniteElementMaterial3D : IFiniteElementMaterial
    {
        double[] Stresses { get; }
        IMatrix2D ConstitutiveMatrix { get; }
        void UpdateMaterial(double[] strains);
        void ClearState();
        void SaveState();
        void ClearStresses();
    }
}
