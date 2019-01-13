using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace ISAAR.MSolve.XFEM.Output.VTK
{
    /// <summary>
    /// Same elements as the original mesh, but each node is included once per element owning it. 
    /// This way we can have multiple field values for the same node.
    /// </summary>
    class VTKMesh2D
    {
        public IReadOnlyList<VTKCell> Cells { get; }
        public IReadOnlyList<VtkPoint2D> Points { get; }

        public VTKMesh2D(IReadOnlyList<VtkPoint2D> points, IReadOnlyList<VTKCell> cells)
        {
            this.Points = points;
            this.Cells = cells;
        }
    }
}
