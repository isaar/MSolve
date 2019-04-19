using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Output
{
    interface IOutputField
    {
        Tensor2D EvaluateAt(XContinuumElement2D element, NaturalPoint point,
            Vector standardDisplacements, Vector enrichedDisplacements);
    }
}
