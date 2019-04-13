using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.Geometry.Tensors;

namespace ISAAR.MSolve.XFEM.Output
{
    interface IOutputField
    {
        Tensor2D EvaluateAt(XContinuumElement2D element, NaturalPoint2D point,
            Vector standardDisplacements, Vector enrichedDisplacements);
    }
}
