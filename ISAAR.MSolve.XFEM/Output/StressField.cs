using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Tensors;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.XFEM.LinearAlgebra;

namespace ISAAR.MSolve.XFEM.Output
{
    class StressField: IOutputField
    {
        public Tensor2D EvaluateAt(XContinuumElement2D element, INaturalPoint2D point, 
            double[] standardDisplacements, double[] enrichedDisplacements)
        {
            EvaluatedInterpolation2D evaluatedInterpolation =
                element.Interpolation.EvaluateAt(element.Nodes, point);
            DenseMatrix displacementGradient = element.CalculateDisplacementFieldGradient(
                point, evaluatedInterpolation, standardDisplacements, enrichedDisplacements);
            Matrix2D constitutive =
                element.Material.CalculateConstitutiveMatrixAt(point, evaluatedInterpolation);

            return element.CalculateStressTensor(displacementGradient, constitutive);
        }
    }
}
