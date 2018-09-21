
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface ICohesiveZoneMaterial3D  
    {
        double[] Tractions { get; }
        IMatrix2D ConstitutiveMatrix { get; }
        void UpdateMaterial(double[] strains);
        int ID { get; }
        bool Modified { get; }
        void ResetModified();
        void SaveState();
        void ClearState();
        void ClearTractions();
        ICohesiveZoneMaterial3D Clone();
    }
}
