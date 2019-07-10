using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Entities.Decomposition
{
    class PolygonalRegion: IRegion2D
    {
        private readonly ConvexPolygon2D polygon;
        private readonly HashSet<LineSegment2D> boundaries;
        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="allVertices"></param>
        /// <param name="boundaryVertices">If the all vertices are boundary, then this list needs to start and end with the same 
        ///     vertex.</param>
        public PolygonalRegion(IReadOnlyList<CartesianPoint> allVertices, HashSet<LineSegment2D> boundaries)
        {
            this.polygon = ConvexPolygon2D.CreateUnsafe(allVertices);
            this.boundaries = boundaries;
        }

        public IReadOnlyList<CartesianPoint> Outline { get { return polygon.Vertices; } }

        public NodePosition FindRelativePosition(XNode node)
        {
            PolygonPointPosition pos = polygon.FindRelativePositionOfPoint(node);
            if (pos == PolygonPointPosition.Outside) return NodePosition.External;
            else if ((pos == PolygonPointPosition.OnEdge) || (pos == PolygonPointPosition.OnVertex))
            {
                // On the inter-subdomain part of the polygon's outline
                foreach (var boundary in boundaries)
                {
                    if (boundary.FindRelativePositionOfPoint(node) == LineSegment2D.SegmentPointPosition.PointOnSegment)
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
