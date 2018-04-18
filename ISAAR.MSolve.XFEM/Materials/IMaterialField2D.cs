using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;

namespace ISAAR.MSolve.XFEM.Materials
{
    // TODO: Perhaps I should return a MaterialPoint DTO instead of each property separately. The trick would be to 
    // avoid redundant object creation and copying doubles in the concrete classes.
    interface IMaterialField2D
    {
        double GetYoungModulusAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation);
        double GetEquivalentYoungModulusAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation);
        double GetPoissonRatioAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation);
        double GetEquivalentPoissonRatioAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation);
        double GetThicknessAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation);
        Matrix CalculateConstitutiveMatrixAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation);

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
