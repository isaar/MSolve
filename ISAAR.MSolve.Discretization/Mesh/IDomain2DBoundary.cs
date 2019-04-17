using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.Discretization.Mesh
{
    public interface IDomain2DBoundary
    {
        // Not on the boundary exactly.
        bool IsInside(CartesianPoint point);
    }
}
