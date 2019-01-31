using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Shapes;

//TODO: this should be moved to Geometry.Shapes once Node2D has been replaced with a generic IPoint2D
namespace ISAAR.MSolve.Preprocessor.Meshes
{
    /// <summary>
    /// Data Transfer Object that packs the <see cref="Geometry.Shapes.CellType"/> with the vertices of a cell. Since there are no 
    /// dependencies, it can be used to transfer cell/element geometry data from one module to another.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CellConnectivity2D
    {
        public CellConnectivity2D(CellType cellType, IReadOnlyList<Node2D> vertices)
        {
            this.CellType = cellType;
            this.Vertices = vertices;
        }

        public CellType CellType { get; }
        public IReadOnlyList<Node2D> Vertices { get; }
    }
}
