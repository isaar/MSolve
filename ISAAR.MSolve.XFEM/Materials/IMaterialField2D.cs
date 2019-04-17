using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Interpolation;

namespace ISAAR.MSolve.XFEM.Materials
{
    // TODO: Perhaps I should return a MaterialPoint DTO instead of each property separately. The trick would be to 
    // avoid redundant object creation and copying doubles in the concrete classes.
    interface IMaterialField2D
    {
        double GetYoungModulusAt(NaturalPoint point, EvaluatedInterpolation2D interpolation);
        double GetEquivalentYoungModulusAt(NaturalPoint point, EvaluatedInterpolation2D interpolation);
        double GetPoissonRatioAt(NaturalPoint point, EvaluatedInterpolation2D interpolation);
        double GetEquivalentPoissonRatioAt(NaturalPoint point, EvaluatedInterpolation2D interpolation);
        double GetThicknessAt(NaturalPoint point, EvaluatedInterpolation2D interpolation);
        Matrix CalculateConstitutiveMatrixAt(NaturalPoint point, EvaluatedInterpolation2D interpolation);

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
