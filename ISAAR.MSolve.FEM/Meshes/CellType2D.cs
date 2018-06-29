using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Meshes
{
    /// <summary>
    /// Defines the shape of a cell only. Since there are no dependencies, it can also be used to map corresponding cell/element 
    /// types from one module to another.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public enum CellType2D
    {
        Quad4, Quad8, Quad9, Tri3, Tri6
    }
}
