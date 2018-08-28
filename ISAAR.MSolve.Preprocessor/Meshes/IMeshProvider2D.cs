using System;
using System.Collections.Generic;
using System.Text;

//TODO: Once FEM.Entities.Node is purged, the generic parameter TVertex should be constrained by IPoint2D. 
//      Also Mesh generation should be moved to Geometry project.
//TODO: Either return a dedicated mesh class (but I need other functionality from meshes...) or make CreateMesh() void and
//      and acces vertices and cells through properties.
namespace ISAAR.MSolve.Preprocessor.Meshes
{
    /// <summary>
    /// Creates 2D meshes for use in FEM or similar methods.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TVertex"></typeparam>
    public interface IMeshProvider2D<TVertex> //where TVertex:IPoint2D 
    {
        (IReadOnlyList<TVertex> vertices, IReadOnlyList<CellConnectivity2D> cells) CreateMesh();
    }
}
