using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh
{
    class CellPlaceholder2D: ICell
    {
        public IReadOnlyList<CartesianPoint2D> Vertices { get; }

        public CellPlaceholder2D(IReadOnlyList<CartesianPoint2D> vertices)
        {
            this.Vertices = vertices;
        }
    }
}
