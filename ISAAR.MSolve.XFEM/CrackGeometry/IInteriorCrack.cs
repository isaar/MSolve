using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Interpolation;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    interface IInteriorCrack
    {
        ICartesianPoint2D StartTip { get; }
        ICartesianPoint2D EndTip { get; }
        List<XContinuumElement2D> StartTipElements { get; }
        List<XContinuumElement2D> EndTipElements { get; }
        TipCoordinateSystem GetTipSystem(CrackTipPosition tip);

        void InitializeGeometry(ICartesianPoint2D startTip, ICartesianPoint2D endTip);
        void UpdateGeometry(double localGrowthAngleStart, double growthLengthStart,
            double localGrowthAngleEnd, double growthLengthEnd); // Perhaps the global angle should be passed in

        double SignedDistanceOf(XNode2D node);
        double SignedDistanceOf(INaturalPoint2D point, XContinuumElement2D element,
            EvaluatedInterpolation2D interpolation);

        IReadOnlyList<TriangleCartesian2D> TriangulateAreaOf(XContinuumElement2D element);
        void UpdateEnrichments();
    }
}
