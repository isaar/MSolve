using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IShellMaterial_v2 : IFiniteElementMaterial
    {
        double[] Stresses { get; }
        IMatrixView ConstitutiveMatrix { get; }
        new IShellMaterial_v2 Clone();
        double[] NormalVectorV3 { get; set; } 
        double[] TangentVectorV1 { get; set; }
        void UpdateMaterial(double[] greenLagrangeStrains);
    }
}
