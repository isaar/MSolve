using System;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Geometry.Shapes
{
    /// <summary>
    /// I would name it Vector, but that means a gazillion of things.
    /// TODO: I don't see why I need the LineSegment2D inheritance; I override pretty much everything.
    /// </summary>
    public class DirectedSegment2D
    {
        /// <summary>
        /// a is the counter-clockwise angle from the global x axis to the local x axis
        /// </summary>
        private readonly double cosa, sina;

        /// <summary>
        /// The coordinates of the global system's origin in the local system
        /// </summary>
        private readonly double originLocalX, originLocalY;

        /// <summary>
        /// The unit vector that is perpendicular to the segment and faces towards the positive local y axis. 
        /// It is constant for a linear segment, so caching it avoids recalculations.
        /// </summary>
        private readonly Vector2 normalVector;

        public CartesianPoint2D Start { get; }
        public CartesianPoint2D End { get; }
        public double Length { get; }
        
        public DirectedSegment2D(CartesianPoint2D start, CartesianPoint2D end)
        {
            this.Start = start;
            this.End = end;

            double startX = start.X;
            double startY = start.Y;
            double dx = end.X - startX;
            double dy = end.Y - startY;

            Length = Math.Sqrt(dx * dx + dy * dy);
            cosa = dx / Length;
            sina = dy / Length;

            originLocalX = -cosa * startX - sina * startY;
            originLocalY = sina * startX - cosa * startY;

            normalVector = Vector2.Create(-sina, cosa ); // This is the opposite from the one in LineSegment2D I think
        }

        public double SignedDistanceOf(CartesianPoint2D point)
        {
            return LocalYOf(point); // This suffices since there is no scaling involved in the transformation.
        }

        public CartesianPoint2D TransformGlobalToLocalPoint(CartesianPoint2D point)
        {
            return new CartesianPoint2D(cosa * point.X + sina * point.Y + originLocalX,
                -sina * point.X + cosa * point.Y + originLocalY);
        }

        // The normal vector for the positive region.
        public Vector2 NormalVectorThrough(CartesianPoint2D point)
        {
            return normalVector.Copy();
        }

        //Perhaps I can override the intersects method and project the other segment onto the local system!
        public LineSegment2D.SegmentSegmentPosition IntersectionWith(LineSegment2D segment, 
            out CartesianPoint2D intersectionPoint)
        {
            return segment.IntersectionWith(new LineSegment2D(Start, End), out intersectionPoint);
        }

        public PointProjectionPosition FindPositionOfProjectionOfPointOntoThis(CartesianPoint2D point)
        {
            double projectionX = LocalXOf(point);
            if (projectionX < 0.0) return PointProjectionPosition.BeforeSegment;
            else if (projectionX > Length) return PointProjectionPosition.AfterSegment;
            else return PointProjectionPosition.OnSegment;
        }

        private double LocalXOf(CartesianPoint2D point)
        {
            return cosa * point.X + sina * point.Y + originLocalX;
        }

        private double LocalYOf(CartesianPoint2D point)
        {
            return -sina * point.X + cosa * point.Y + originLocalY;
        }

        public enum PointProjectionPosition
        {
            BeforeSegment, OnSegment, AfterSegment
        }
    }
}
