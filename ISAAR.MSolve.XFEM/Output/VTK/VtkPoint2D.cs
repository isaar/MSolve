using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    internal class VtkPoint2D: CartesianPoint2D
    {
        public VtkPoint2D(int id, CartesianPoint2D originalNode): base(originalNode.X, originalNode.Y)
        {
            this.ID = id;
        }

        public int ID { get; }
    }
}
