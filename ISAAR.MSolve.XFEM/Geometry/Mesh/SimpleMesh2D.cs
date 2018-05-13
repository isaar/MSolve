using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh
{
    /// <summary>
    /// Unoptimized. All search operations take linear time: O(nodesCount) and O(elementsCount).
    /// </summary>
    /// <typeparam name="TVertex"></typeparam>
    /// <typeparam name="TCell"></typeparam>
    class SimpleMesh2D<TVertex, TCell>: IMesh2D<TVertex, TCell> 
        where TVertex: ICartesianPoint2D 
        where TCell: class, ICell
    {
        private readonly IDomainBoundary boundary;

        public IReadOnlyList<TVertex> Vertices { get; }
        public IReadOnlyList<TCell> Cells { get; }

        public SimpleMesh2D(IReadOnlyList<TVertex> vertices, IReadOnlyList<TCell> faces, IDomainBoundary boundary)
        {
            this.Vertices = vertices;
            this.Cells = faces;
            this.boundary = boundary;
        }

        // TODO: handle cases where more the point lies on an element edge or node.
        public IReadOnlyList<TCell> FindElementsContainingPoint(ICartesianPoint2D point, TCell startingElement = null)
        {
            var containingElements = new List<TCell>();
            foreach (TCell element in Cells) // O(elementsCount)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Vertices);
                PolygonPointPosition pos = outline.FindRelativePositionOfPoint(point);
                if ((pos == PolygonPointPosition.Inside) || (pos == PolygonPointPosition.OnEdge) ||
                    (pos == PolygonPointPosition.OnVertex)) containingElements.Add(element);
            }
            return containingElements;
        }

        public IReadOnlyList<TCell> FindElementsIntersectedByCircle(Circle2D circle, TCell startingElement = null)
        {
            var intersectedElements = new List<TCell>();
            foreach (TCell element in Cells)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Vertices);
                CirclePolygonPosition relativePosition = outline.FindRelativePositionOfCircle(circle);
                if (relativePosition == CirclePolygonPosition.Intersecting) intersectedElements.Add(element);
            }
            return intersectedElements;
        }

        public IReadOnlyList<TCell> FindElementsInsideCircle(Circle2D circle, TCell startingElement = null)
        {
            var internalElements = new List<TCell>();
            foreach (TCell element in Cells)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Vertices);
                CirclePolygonPosition relativePosition = outline.FindRelativePositionOfCircle(circle);
                if ((relativePosition == CirclePolygonPosition.Intersecting) ||
                    (relativePosition == CirclePolygonPosition.PolygonInsideCircle)) internalElements.Add(element);
            }
            return internalElements;
        }

        public ISet<TCell> FindElementsWithNode(TVertex node)
        {
            var neighboringElements = new HashSet<TCell>();
            foreach (TCell element in Cells)
            {
                if (element.Vertices.Contains(node)) neighboringElements.Add(element);
            }
            return neighboringElements;
        }

        public IReadOnlyList<TVertex> FindNodesInsideCircle(Circle2D circle,
            bool findBoundaryNodes = true, TCell startingElement = null)
        {
            var selectedNodes = new List<TVertex>();
            foreach (TVertex node in Vertices) // O(nodesCount)
            {
                CirclePointPosition relativePosition = circle.FindRelativePositionOfPoint(node);
                if ( (relativePosition == CirclePointPosition.Inside) ||
                    (findBoundaryNodes && (relativePosition == CirclePointPosition.On)) )
                {
                    selectedNodes.Add(node);
                }
            }
            return selectedNodes;
        }

        public bool IsInsideBoundary(ICartesianPoint2D point)
        {
            return boundary.IsInside(point);
        }
    }
}
