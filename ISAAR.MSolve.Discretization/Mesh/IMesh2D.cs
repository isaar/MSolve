using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.Discretization.Mesh
{
    public interface IMesh2D<TNode, TElement> 
        where TNode: INode 
        where TElement: class, ICell<TNode>
    {
        IReadOnlyList<TNode> Nodes { get; }
        IReadOnlyList<TElement> Elements { get; }

        /// <summary>
        /// Find the elements that contain the provided point. The following cases are possible:
        /// Case 1) The point is inside 1 element. Only that element will be returned. 
        /// Case 2) The point is on the edge of 2 elements. Both elements will be returned.
        /// Case 3) The point coincides with a node of 3 or more elements. All these elements will be returned.
        /// Case 4) The point is outside the mesh. The returned collection will be empty.
        /// </summary>
        /// <param name="point">The point to search.</param>
        /// <param name="startingElement">An element around which to start the search. If no such element is provided 
        ///     the search will probably take longer, usually O(elementsCount).</param>
        /// <returns></returns>
        IReadOnlyList<TElement> FindElementsContainingPoint(CartesianPoint point, TElement startingElement = null);

        /// <summary>
        /// Find the elements that are intersected by the provided circle. The following cases are possible:
        /// Case 1) The circle intersects the edges of some elements and at least 1 intesection point per element is 
        ///     not a node. These elements will be returned.
        /// Case 2) All elements are either completely outside or completely inside the circle (they might touch at 
        ///     nodal points though). A 3rd option is that the circle is small enough to be completely inside an 
        ///     element. The returned collection will be empty.
        /// </summary>
        /// <param name="circle">The circle to consider.</param>
        /// <param name="startingElement">An element around which to start the search. If no such element is provided 
        ///     the search will probably take longer, usually O(elementsCount).</param>
        /// <returns></returns>
        IReadOnlyList<TElement> FindElementsIntersectedByCircle(Circle2D circle, TElement startingElement = null);

        /// <summary>
        /// Find the elements that are completely inside or intersected by the provided circle. 
        /// <param name="circle">The circle to consider.</param>
        /// <param name="startingElement">An element around which to start the search. If no such element is provided 
        ///     the search will probably take longer, usually O(elementsCount).</param>
        /// <returns></returns>
        IReadOnlyList<TElement> FindElementsInsideCircle(Circle2D circle, TElement startingElement = null);

        ISet<TElement> FindElementsWithNode(TNode node);

        /// <summary>
        /// Find the nodes of the mesh which are inside (and optionally exactly on) the provided circle. If there are
        /// none, the returned collection will be empty.
        /// </summary>
        /// <param name="circle"></param>
        /// <param name="findBoundaryNodes">True if nodes, whose distance to the circle's center is exactly equal to 
        ///     the radius, are to be returned as well.</param>
        /// <param name="startingElement">An element around which to start the search. If no such element is provided 
        ///     the search will probably take longer, usually O(nodesCount).</param>
        /// <returns></returns>
        IReadOnlyList<TNode> FindNodesInsideCircle(Circle2D circle, 
            bool findBoundaryNodes = true, TElement startingElement = null);

        bool IsInsideBoundary(CartesianPoint point);
    }
}
