using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Output
{
    interface IOutputField
    {
        Tensor2D EvaluateAt(XContinuumElement2D element, INaturalPoint2D point,
            double[] standardDisplacements, double[] enrichedDisplacements);
    }
}
