using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Geometry.Boundaries
{
    class RectangularBoundary: IDomainBoundary
    {
        private readonly double minX;
        private readonly double maxX;
        private readonly double minY;
        private readonly double maxY;

        public RectangularBoundary(double minX, double maxX, double minY, double maxY)
        {
            this.minX = minX;
            this.maxX = maxX;
            this.minY = minY;
            this.maxY = maxY;
        }

        public bool IsInside(CartesianPoint2D point)
        {
            if ((point.X > minX) && (point.X < maxX) && (point.Y > minY) && (point.Y < maxY)) return true;
            else return false;
        }
    }
}
