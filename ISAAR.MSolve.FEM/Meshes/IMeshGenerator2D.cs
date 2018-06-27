using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;

//TODO: Once FEM.Entities.Node is purged, the generic parameter TVertex should be constrained by IPoint2D. 
//      Also Mesh generation should be moved to Geometry project.
namespace ISAAR.MSolve.FEM.Meshes
{
    /// <summary>
    /// Creates 2D meshes for use in FEM or similar methods.
    /// </summary>
    /// <typeparam name="TVertex"></typeparam>
    public interface IMeshGenerator2D<TVertex> //where TVertex:IPoint2D 
    {
        (TVertex[] vertices, TVertex[][] cellConnectivity) CreateMesh();
    }
}
