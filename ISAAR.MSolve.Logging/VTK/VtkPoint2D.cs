using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging.VTK
{
    /// <summary>
    /// Vertex used to represent VTK grids.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class VtkPoint2D: CartesianPoint2D
    {
        public VtkPoint2D(int id, double x, double y): base(x, y)
        {
            this.ID = id;
        }

        public int ID { get; }
    }
}
