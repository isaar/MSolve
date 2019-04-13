using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Entities.Decomposition
{
    class RectangularRegion: IRegion2D
    {
        private readonly double minX;
        private readonly double extendedMinX;
        private readonly double maxX;
        private readonly double extendedMaxX;
        private readonly double minY;
        private readonly double extendedMinY;
        private readonly double maxY;
        private readonly double extendedMaxY;
        private readonly double tolerance;
        private readonly Dictionary<RectangleEdge, bool> boundaries;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="minX"></param>
        /// <param name="maxX"></param>
        /// <param name="minY"></param>
        /// <param name="maxY"></param>
        /// <param name="tolerance">How close to the boundary a node must be to consider it as a boundary node. If you set it to 
        ///     0.0, it will check whether the corrdinates are exactly the same.</param>
        public RectangularRegion(double minX, double maxX, double minY, double maxY, double tolerance)
        {
            this.minX = minX;
            this.extendedMinX = minX - tolerance;
            this.maxX = maxX;
            this.extendedMaxX = maxX + tolerance;
            this.minY = minY;
            this.extendedMinY = minY - tolerance;
            this.maxY = maxY;
            this.extendedMaxY = maxY + tolerance;
            this.tolerance = tolerance;

            boundaries = new Dictionary<RectangleEdge, bool>();
            boundaries.Add(RectangleEdge.Left, false);
            boundaries.Add(RectangleEdge.Right, false);
            boundaries.Add(RectangleEdge.Up, false);
            boundaries.Add(RectangleEdge.Down, false);
        }

        public void AddBoundaryEdge(RectangleEdge edge)
        {
            boundaries[edge] = true;
        }

        public NodePosition FindRelativePosition(XNode2D node)
        {
            //TODO: there might be a faster way to check these
            if ( (node.X >= extendedMinX) && (node.X <= extendedMaxX) )
            {
                if ((node.Y >= extendedMinY) && (node.Y <= extendedMaxY))
                {
                    if (IsBoundaryLeft(node) || IsBoundaryRight(node) || IsBoundaryUp(node) || IsBoundaryDown(node))
                    { //TODO: The OOP way would be to have an enum class that answers, but that is far more complex here.
                        return NodePosition.Boundary; 
                    }
                    else return NodePosition.Internal;
                }
            }
            return NodePosition.External;
        }

        private bool IsBoundaryLeft(XNode2D node)
        {
            if (boundaries[RectangleEdge.Left]) return Math.Abs(node.X - minX) <= tolerance;
            return false;
        }

        private bool IsBoundaryRight(XNode2D node)
        {
            if (boundaries[RectangleEdge.Right]) return Math.Abs(node.X - maxX) <= tolerance;
            return false;
        }

        private bool IsBoundaryUp(XNode2D node)
        {
            if (boundaries[RectangleEdge.Up]) return Math.Abs(node.Y - maxY) <= tolerance;
            return false;
        }

        private bool IsBoundaryDown(XNode2D node)
        {
            if (boundaries[RectangleEdge.Down]) return Math.Abs(node.Y - minY) <= tolerance;
            return false;
        }

        public enum RectangleEdge
        {
            Left, Right, Up, Down
        }
    }
}
