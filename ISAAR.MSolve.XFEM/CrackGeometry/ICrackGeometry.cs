using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Interpolation;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    interface ICrackGeometry
    {
        double SignedDistanceOf(XNode2D node);
        double SignedDistanceOf(INaturalPoint2D point, XContinuumElement2D element,
            EvaluatedInterpolation2D interpolation);

        ICartesianPoint2D GetCrackTip(CrackTipPosition tipPosition);
        /// TODO: an interface is needed for TipSystems. Then the explicit (global, local, polar) systems or the level 
        /// sets could be used for the transformations (points, vectors, derivatives)
        TipCoordinateSystem GetTipSystem(CrackTipPosition tipPosition);
        IReadOnlyList<XContinuumElement2D> GetTipElements(CrackTipPosition tipPosition);

        SortedSet<ICartesianPoint2D> FindTriangleVertices(XContinuumElement2D element);
        void UpdateEnrichments();
    }
}
