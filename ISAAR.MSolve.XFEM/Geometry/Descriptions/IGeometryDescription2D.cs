using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Descriptions
{
    interface IGeometryDescription2D
    {
        double SignedDistanceOf(ICartesianPoint2D point);
        Vector2 NormalVectorThrough(ICartesianPoint2D point);

        ICartesianPoint2D StartPoint { get; }

        // This one's orientation requires more thought, especially since the convention for determining the
        // level set value gives the opposite sign from the rest of the curve
        //double StartPointOrientation() { get; } 

        ICartesianPoint2D EndPoint { get; }

        /// <summary>
        /// Counter-clockwise angle from global cartesian x axis to a vector which 1) starts at the end point of the 
        /// curve, 2) is tangent to the curve and 3) heads outwards from the curve.
        /// TODO: perhaps it should return a local coordinate system. I do not need the angle, but its cos and sign.
        /// </summary>
        double EndPointOrientation();

        double StartPointOrientation();

        // Perhaps geometry classes should be decoupled from elements and interact through polygons instead.
        IReadOnlyList<ICartesianPoint2D> IntersectionWith(XContinuumElement2D element);

        // This is the correct one but it needs constrained Delauny triangulation
        //public void IntersectionWith(ContinuumElement2D element, out ICartesianPoint2D[] points, out LineSegment2D[] segments);
    }
}
