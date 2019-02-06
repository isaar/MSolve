using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IStochasticIsotropicContinuumMaterial3D_v2: IIsotropicContinuumMaterial3D_v2 //maybe we don't need it?? the interface members are unused
    {
        IStochasticMaterialCoefficientsProvider CoefficientsProvider { get; set; }
        IMatrixView GetConstitutiveMatrix(double[] coordinates);
    }
}
