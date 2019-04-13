using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Geometry.Shapes
{
    public interface IOpenCurve2D : ICurve2D
    {
        CartesianPoint2D Start { get; }
        CartesianPoint2D End { get; }

        /// <summary>
        /// Unit vector. It will coincide with the normal vector if rotated -PI/2.
        /// </summary>
        Vector2 TangentAtStart { get; }

        /// <summary>
        /// Unit vector. It will coincide with the normal vector if rotated PI/2.
        /// </summary>
        Vector2 TangentAtEnd { get; }
    }
}
