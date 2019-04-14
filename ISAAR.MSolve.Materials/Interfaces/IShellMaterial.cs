using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IShellMaterial : IFiniteElementMaterial
    {
        double[] Stresses { get; }
        IMatrixView ConstitutiveMatrix { get; }
        new IShellMaterial Clone();
        double[] NormalVectorV3 { get; set; }
        double[] TangentVectorV2 { get; set; }
        double[] TangentVectorV1 { get; set; }
        void UpdateMaterial(double[] greenLagrangeStrains);
    }
}
