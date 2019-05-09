using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

//TODO: make this class generic on TVertex and TCell. To do so, TCell must also be generic on TVertex, which creates a whole
//      chain of difficulties.
namespace ISAAR.MSolve.Discretization.Mesh
{
    /// <summary>
    /// Also keeps track of which elements are connected to each node
    /// </summary>
    public class BidirectionalMesh2D<TNode, TElement>: IMesh2D<TNode, TElement>
        where TNode : INode
        where TElement : class, ICell<TNode>
    {
        private readonly IDomain2DBoundary boundary;
        private readonly IReadOnlyDictionary<TNode, HashSet<TElement>> nodeConnectivity;

        public IReadOnlyList<TNode> Nodes { get; }
        public IReadOnlyList<TElement> Elements { get; }

        public BidirectionalMesh2D(IReadOnlyList<TNode> vertices, IReadOnlyList<TElement> faces, IDomain2DBoundary boundary)
        {
            this.Nodes = vertices;
            this.Elements = faces;
            this.boundary = boundary;

            // Connect cells to vertices
            var connectivity = new Dictionary<TNode, HashSet<TElement>>();
            foreach (var cell in Elements)
            {
                foreach (var vertex in cell.Nodes)
                {
                    bool exists = connectivity.TryGetValue(vertex, out HashSet<TElement> vertexCells);
                    if (exists) vertexCells.Add(cell);
                    else
                    {
                        vertexCells = new HashSet<TElement>();
                        vertexCells.Add(cell);
                        connectivity.Add(vertex, vertexCells);
                    }
                }
            }
            this.nodeConnectivity = connectivity;
        }

        // TODO: handle cases where more the point lies on an element edge or node.
        // TODO: Do it by spreading around the starting cell, breadth-first-search
        public IReadOnlyList<TElement> FindElementsContainingPoint(CartesianPoint point,
            TElement startingElement = null)
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

        // TODO: Do it by spreading around the starting cell, breadth-first-search
        public IReadOnlyList<TElement> FindElementsIntersectedByCircle(Circle2D circle,
            TElement startingElement = null)
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

        // TODO: Do it by spreading around the starting cell, breadth-first-search
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
            return new HashSet<TElement>(nodeConnectivity[node]);
        }

        // TODO: Do it by spreading around the starting cell, breadth-first-search
        public IReadOnlyList<TNode> FindNodesInsideCircle(Circle2D circle, 
            bool findBoundaryNodes = true, TElement startingElement = null)
        {
            var selectedNodes = new List<TNode>();
            foreach (var node in Nodes) // O(nodesCount)
            {
                CirclePointPosition relativePosition = circle.FindRelativePositionOfPoint(new CartesianPoint(node.X, node.Y));
                if ((relativePosition == CirclePointPosition.Inside) ||
                    (findBoundaryNodes && (relativePosition == CirclePointPosition.On)))
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
