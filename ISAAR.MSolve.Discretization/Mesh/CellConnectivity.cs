using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;

//TODO: this should be moved to Geometry.Shapes once Node2D has been replaced with a generic IPoint2D
namespace ISAAR.MSolve.Discretization.Mesh
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
