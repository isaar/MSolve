using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Meshes
{
    /// <summary>
    /// Data Transfer Object that packs the <see cref="CellType2D"/> with the vertices of a cell. Since there are no 
    /// dependencies, it can be used to transfer cell/element geometry data from one module to another.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CellConnectivity2D
    {
        public CellConnectivity2D(CellType2D cellType, IReadOnlyList<Node2D> vertices)
        {
            this.CellType = cellType;
            this.Vertices = vertices;
        }

        public CellType2D CellType { get; }
        public IReadOnlyList<Node2D> Vertices { get; }
    }
}
