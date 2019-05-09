using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Geometry.Shapes
{
    public interface ICurve2D
    {
        double SignedDistanceOf(CartesianPoint point);

        /// <summary>
        /// Unit length
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        Vector2 NormalVectorThrough(CartesianPoint point);
    }
}
