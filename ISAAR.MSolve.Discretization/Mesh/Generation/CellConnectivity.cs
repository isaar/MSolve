using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: this should be moved to Geometry.Shapes once Node2D has been replaced with a generic IPoint
namespace ISAAR.MSolve.Discretization.Mesh.Generation
{
    /// <summary>
    /// Data Transfer Object that packs the <see cref="CellType"/> with the vertices of a cell. Since there are no 
    /// dependencies, it can be used to transfer cell/element geometry data from one module to another.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CellConnectivity<TNode> where TNode: INode
    {
        public CellConnectivity(CellType cellType, IReadOnlyList<TNode> vertices)
        {
            this.CellType = cellType;
            this.Vertices = vertices;
        }

        public CellType CellType { get; }
        public IReadOnlyList<TNode> Vertices { get; }
    }
}
