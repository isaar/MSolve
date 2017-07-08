using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh
{
    class MeshFacePlaceholder2D: IMeshFace
    {
        public IReadOnlyList<ICartesianPoint2D> Vertices { get; }

        public MeshFacePlaceholder2D(IReadOnlyList<ICartesianPoint2D> vertices)
        {
            this.Vertices = vertices;
        }
    }
}
