using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.CrackPropagation
{
    interface IPropagator
    {
        PropagationLogger Logger { get; }

        (double growthAngle, double growthLength) Propagate(
            IDofOrderer dofOrderer, Vector totalFreeDisplacements, Vector totalConstrainedDisplacements,
            ICartesianPoint2D crackTip, TipCoordinateSystem tipSystem, IReadOnlyList<XContinuumElement2D> tipElements);
    }
}
