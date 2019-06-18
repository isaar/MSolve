using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning
{
    public class PolygonalRegion2D : IRegion2D
    {
        private readonly ConvexPolygon2D polygon;
        private readonly IEnumerable<LineSegment2D> boundaries;
        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="allVertices"></param>
        /// <param name="boundaryVertices">If the all vertices are boundary, then this list needs to start and end with the same 
        ///     vertex.</param>
        public PolygonalRegion2D(IReadOnlyList<CartesianPoint> allVertices, IEnumerable<LineSegment2D> boundaries)
        {
            this.polygon = ConvexPolygon2D.CreateUnsafe(allVertices);
            this.boundaries = boundaries;
        }

        public IReadOnlyList<CartesianPoint> Outline { get { return polygon.Vertices; } }

        public NodePosition FindRelativePosition(INode node)
        {
            var point = new CartesianPoint(node.X, node.Y, node.Z);
            PolygonPointPosition pos = polygon.FindRelativePositionOfPoint(point);
            if (pos == PolygonPointPosition.Outside) return NodePosition.External;
            else if ((pos == PolygonPointPosition.OnEdge) || (pos == PolygonPointPosition.OnVertex))
            {
                // On the inter-subdomain part of the polygon's outline
                foreach (var boundary in boundaries)
                {
                    if (boundary.FindRelativePositionOfPoint(point) == LineSegment2D.SegmentPointPosition.PointOnSegment)
                    {
                        return NodePosition.Boundary;
                    }
                }

                // otherwise
                return NodePosition.Internal;
            }
            else if (pos==PolygonPointPosition.Inside)
            {
                return NodePosition.Internal;
            }
            else
            {
                throw new Exception("This code should not have been reached");
            }
        }
    }
}
