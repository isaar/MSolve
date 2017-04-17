using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh
{
    /// <summary>
    /// Unoptimized. All search operations take linear time: O(nodesCount) and O(elementsCount).
    /// </summary>
    /// <typeparam name="TVertex"></typeparam>
    /// <typeparam name="TFace"></typeparam>
    class SimpleMesh2D<TVertex, TFace>: IMesh2D<TVertex, TFace> 
        where TVertex: ICartesianPoint2D 
        where TFace: class, IMeshFace
    {
        public IReadOnlyList<TVertex> Vertices { get; }
        public IReadOnlyList<TFace> Faces { get; }

        public SimpleMesh2D(IReadOnlyList<TVertex> vertices, IReadOnlyList<TFace> faces)
        {
            this.Vertices = vertices;
            this.Faces = faces;
        }

        // TODO: handle cases where more the point lies on an element edge or node.
        public IReadOnlyList<TFace> FindElementsContainingPoint(ICartesianPoint2D point, TFace startingElement = null)
        {
            var containingElements = new List<TFace>();
            foreach (TFace element in Faces) // O(elementsCount)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Vertices);
                PolygonPointPosition pos = outline.FindRelativePositionOfPoint(point);
                if ((pos == PolygonPointPosition.Inside) || (pos == PolygonPointPosition.OnEdge) ||
                    (pos == PolygonPointPosition.OnVertex)) containingElements.Add(element);
            }
            return containingElements;
        }

        public IReadOnlyList<TFace> FindElementsIntersectedByCircle(Circle2D circle, TFace startingElement = null)
        {
            var intersectedElements = new List<TFace>();
            foreach (TFace element in Faces)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Vertices);
                CirclePolygonPosition relativePosition = outline.FindRelativePositionOfCircle(circle);
                if (relativePosition == CirclePolygonPosition.Intersecting) intersectedElements.Add(element);
            }
            return intersectedElements;
        }

        public IReadOnlyList<TVertex> FindNodesInsideCircle(Circle2D circle,
            bool findBoundaryNodes = true, TFace startingElement = null)
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

    }
}
