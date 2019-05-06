using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Once FEM.Entities.Node is purged, the generic parameter TVertex should be constrained by IPoint. 
//      Also Mesh generation should be moved to Geometry project.
//TODO: Either return a dedicated mesh class (but I need other functionality from meshes...) or make CreateMesh() void and
//      and acces vertices and cells through properties.
namespace ISAAR.MSolve.Discretization.Mesh.Generation
{
    /// <summary>
    /// Creates 2D and 3D meshes for use in FEM or similar methods.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TNode"></typeparam>
    public interface IMeshGenerator<TNode> where TNode : INode 
    {
        (IReadOnlyList<TNode> nodes, IReadOnlyList<CellConnectivity<TNode>> elements) CreateMesh(CreateNode<TNode> createNode);
    }

    /// <summary>
    /// Initializes a new instance of some class that implements <see cref="INode"/>.
    /// </summary>
    /// <param name="id">The unique identifier of the node.</param>
    /// <param name="x">The x coordinate of the node.</param>
    /// <param name="y">The y coordinate of the node. Leave it 0.0 for 1D problems.</param>
    /// <param name="z">The z coordinate of the node. Leave it 0.0 for 1D, 2D problems.</param>
    public delegate TNode CreateNode<TNode>(int id, double x, double y, double z) where TNode: INode;
}
