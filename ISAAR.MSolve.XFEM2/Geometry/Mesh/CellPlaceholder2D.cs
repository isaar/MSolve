using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh
{
    class CellPlaceholder2D: ICell
    {
        public IReadOnlyList<ICartesianPoint2D> Vertices { get; }

        public CellPlaceholder2D(IReadOnlyList<ICartesianPoint2D> vertices)
        {
            this.Vertices = vertices;
        }
    }
}
