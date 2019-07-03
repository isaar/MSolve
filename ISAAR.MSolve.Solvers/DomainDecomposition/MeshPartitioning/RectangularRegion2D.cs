using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning
{
    public class RectangularRegion2D : IRegion2D
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
        ///     0.0, it will check whether the coordinates are exactly the same.</param>
        public RectangularRegion2D(double minX, double minY, double maxX, double maxY, double tolerance)
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

        public NodePosition FindRelativePosition(INode node)
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

        private bool IsBoundaryLeft(INode node)
        {
            if (boundaries[RectangleEdge.Left]) return Math.Abs(node.X - minX) <= tolerance;
            return false;
        }

        private bool IsBoundaryRight(INode node)
        {
            if (boundaries[RectangleEdge.Right]) return Math.Abs(node.X - maxX) <= tolerance;
            return false;
        }

        private bool IsBoundaryUp(INode node)
        {
            if (boundaries[RectangleEdge.Up]) return Math.Abs(node.Y - maxY) <= tolerance;
            return false;
        }

        private bool IsBoundaryDown(INode node)
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
