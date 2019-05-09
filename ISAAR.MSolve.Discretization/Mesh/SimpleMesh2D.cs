using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.Discretization.Mesh
{
    /// <summary>
    /// Unoptimized. Uses little storage. All search operations take linear time: O(nodesCount) and O(elementsCount).
    /// </summary>
    /// <typeparam name="TNode"></typeparam>
    /// <typeparam name="TElement"></typeparam>
    public class SimpleMesh2D<TNode, TElement>: IMesh2D<TNode, TElement>
        where TNode : INode
        where TElement : class, ICell<TNode>
    {
        private readonly IDomain2DBoundary boundary;

        public IReadOnlyList<TNode> Nodes { get; }
        public IReadOnlyList<TElement> Elements { get; }

        public SimpleMesh2D(IReadOnlyList<TNode> vertices, IReadOnlyList<TElement> faces, IDomain2DBoundary boundary)
        {
            this.Nodes = vertices;
            this.Elements = faces;
            this.boundary = boundary;
        }

        // TODO: handle cases where more the point lies on an element edge or node.
        public IReadOnlyList<TElement> FindElementsContainingPoint(CartesianPoint point, TElement startingElement = null)
        {
            var containingElements = new List<TElement>();
            foreach (TElement element in Elements) // O(elementsCount)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes.ToCartesianPoints());
                PolygonPointPosition pos = outline.FindRelativePositionOfPoint(point);
                if ((pos == PolygonPointPosition.Inside) || (pos == PolygonPointPosition.OnEdge) ||
                    (pos == PolygonPointPosition.OnVertex)) containingElements.Add(element);
            }
            return containingElements;
        }

        public IReadOnlyList<TElement> FindElementsIntersectedByCircle(Circle2D circle, TElement startingElement = null)
        {
            var intersectedElements = new List<TElement>();
            foreach (TElement element in Elements)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes.ToCartesianPoints());
                CirclePolygonPosition relativePosition = outline.FindRelativePositionOfCircle(circle);
                if (relativePosition == CirclePolygonPosition.Intersecting) intersectedElements.Add(element);
            }
            return intersectedElements;
        }

        public IReadOnlyList<TElement> FindElementsInsideCircle(Circle2D circle, TElement startingElement = null)
        {
            var internalElements = new List<TElement>();
            foreach (TElement element in Elements)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes.ToCartesianPoints());
                CirclePolygonPosition relativePosition = outline.FindRelativePositionOfCircle(circle);
                if ((relativePosition == CirclePolygonPosition.Intersecting) ||
                    (relativePosition == CirclePolygonPosition.PolygonInsideCircle)) internalElements.Add(element);
            }
            return internalElements;
        }

        public ISet<TElement> FindElementsWithNode(TNode node)
        {
            var neighboringElements = new HashSet<TElement>();
            foreach (TElement element in Elements)
            {
                if (element.Nodes.Contains(node)) neighboringElements.Add(element);
            }
            return neighboringElements;
        }

        public IReadOnlyList<TNode> FindNodesInsideCircle(Circle2D circle,
            bool findBoundaryNodes = true, TElement startingElement = null)
        {
            var selectedNodes = new List<TNode>();
            foreach (TNode node in Nodes) // O(nodesCount)
            {
                CirclePointPosition relativePosition = circle.FindRelativePositionOfPoint(new CartesianPoint(node.X, node.Y));
                if ( (relativePosition == CirclePointPosition.Inside) ||
                    (findBoundaryNodes && (relativePosition == CirclePointPosition.On)) )
                {
                    selectedNodes.Add(node);
                }
            }
            return selectedNodes;
        }

        public bool IsInsideBoundary(CartesianPoint point)
        {
            return boundary.IsInside(point);
        }
    }
}
