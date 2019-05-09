using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.XFEM.Materials
{
    // TODO: Perhaps I should return a MaterialPoint DTO instead of each property separately. The trick would be to 
    // avoid redundant object creation and copying doubles in the concrete classes.
    public interface IMaterialField2D
    {
        IMaterialField2D Clone();

        double GetYoungModulusAt(NaturalPoint point, EvalInterpolation2D interpolation);
        double GetEquivalentYoungModulusAt(NaturalPoint point, EvalInterpolation2D interpolation);
        double GetPoissonRatioAt(NaturalPoint point, EvalInterpolation2D interpolation);
        double GetEquivalentPoissonRatioAt(NaturalPoint point, EvalInterpolation2D interpolation);
        double GetThicknessAt(NaturalPoint point, EvalInterpolation2D interpolation);
        Matrix CalculateConstitutiveMatrixAt(NaturalPoint point, EvalInterpolation2D interpolation);


        /// <summary>
        /// Used in material stochasticity, XFEM inclusions, etc
        /// </summary>
        void UpdateDistributions();

        /// <summary>
        /// Used in material non-linearity.
        /// </summary>
        void UpdateStrains();
    }
}
