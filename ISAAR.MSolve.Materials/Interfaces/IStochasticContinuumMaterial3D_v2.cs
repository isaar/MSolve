using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IStochasticContinuumMaterial3D_v2 : IContinuumMaterial3D_v2
    {
        IStochasticMaterialCoefficientsProvider CoefficientsProvider { get; set; }
        IMatrixView GetConstitutiveMatrix(double[] coordinates);
    }
}
