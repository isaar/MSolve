using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

//TODO: make this class generic on TVertex and TCell. To do so, TCell must also be generic on TVertex, which creates a whole
//      chain of difficulties.
namespace ISAAR.MSolve.XFEM.Geometry.Mesh
{
    /// <summary>
    /// Also keeps track of which elements are connected to each node
    /// </summary>
    class BiMesh2D : IMesh2D<XNode2D, XContinuumElement2D>
    {
        private readonly IDomainBoundary boundary;
        private readonly IReadOnlyDictionary<XNode2D, HashSet<XContinuumElement2D>> vertexConnectivity;

        public IReadOnlyList<XNode2D> Vertices { get; }
        public IReadOnlyList<XContinuumElement2D> Cells { get; }

        public BiMesh2D(IReadOnlyList<XNode2D> vertices, IReadOnlyList<XContinuumElement2D> faces, IDomainBoundary boundary)
        {
            this.Vertices = vertices;
            this.Cells = faces;
            this.boundary = boundary;

            // Connect cells to vertices
            var connectivity = new Dictionary<XNode2D, HashSet<XContinuumElement2D>>();
            foreach (var cell in Cells)
            {
                foreach (var vertex in cell.Nodes)
                {
                    bool exists = connectivity.TryGetValue(vertex, out HashSet<XContinuumElement2D> vertexCells);
                    if (exists) vertexCells.Add(cell);
                    else
                    {
                        vertexCells = new HashSet<XContinuumElement2D>();
                        vertexCells.Add(cell);
                        connectivity.Add(vertex, vertexCells);
                    }
                }
            }
            this.vertexConnectivity = connectivity;
        }

        // TODO: handle cases where more the point lies on an element edge or node.
        // TODO: Do it by spreading around the starting cell, breadth-first-search
        public IReadOnlyList<XContinuumElement2D> FindElementsContainingPoint(ICartesianPoint2D point, 
            XContinuumElement2D startingElement = null)
        {
            var containingElements = new List<XContinuumElement2D>();
            foreach (XContinuumElement2D element in Cells) // O(elementsCount)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Vertices);
                PolygonPointPosition pos = outline.FindRelativePositionOfPoint(point);
                if ((pos == PolygonPointPosition.Inside) || (pos == PolygonPointPosition.OnEdge) ||
                    (pos == PolygonPointPosition.OnVertex)) containingElements.Add(element);
            }
            return containingElements;
        }

        // TODO: Do it by spreading around the starting cell, breadth-first-search
        public IReadOnlyList<XContinuumElement2D> FindElementsIntersectedByCircle(Circle2D circle, 
            XContinuumElement2D startingElement = null)
        {
            var intersectedElements = new List<XContinuumElement2D>();
            foreach (XContinuumElement2D element in Cells)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Vertices);
                CirclePolygonPosition relativePosition = outline.FindRelativePositionOfCircle(circle);
                if (relativePosition == CirclePolygonPosition.Intersecting) intersectedElements.Add(element);
            }
            return intersectedElements;
        }

        // TODO: Do it by spreading around the starting cell, breadth-first-search
        public IReadOnlyList<XContinuumElement2D> FindElementsInsideCircle(Circle2D circle, 
            XContinuumElement2D startingElement = null)
        {
            var internalElements = new List<XContinuumElement2D>();
            foreach (XContinuumElement2D element in Cells)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Vertices);
                CirclePolygonPosition relativePosition = outline.FindRelativePositionOfCircle(circle);
                if ((relativePosition == CirclePolygonPosition.Intersecting) ||
                    (relativePosition == CirclePolygonPosition.PolygonInsideCircle)) internalElements.Add(element);
            }
            return internalElements;
        }

        public ISet<XContinuumElement2D> FindElementsWithNode(XNode2D node)
        {
            return new HashSet<XContinuumElement2D>(vertexConnectivity[node]);
        }

        // TODO: Do it by spreading around the starting cell, breadth-first-search
        public IReadOnlyList<XNode2D> FindNodesInsideCircle(Circle2D circle,
            bool findBoundaryNodes = true, XContinuumElement2D startingElement = null)
        {
            var selectedNodes = new List<XNode2D>();
            foreach (var node in Vertices) // O(nodesCount)
            {
                CirclePointPosition relativePosition = circle.FindRelativePositionOfPoint(node);
                if ((relativePosition == CirclePointPosition.Inside) ||
                    (findBoundaryNodes && (relativePosition == CirclePointPosition.On)))
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
