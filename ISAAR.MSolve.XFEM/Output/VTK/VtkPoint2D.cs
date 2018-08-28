using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    internal class VtkPoint2D: ICartesianPoint2D
    {
        public VtkPoint2D(int id, ICartesianPoint2D originalNode)
        {
            this.ID = id;
            this.X = originalNode.X;
            this.Y = originalNode.Y;
        }

        public Vector2 Coordinates { get { return Vector2.Create( X, Y ); } }
        public int ID { get; }
        public double X { get; }
        public double Y { get; }
    }
}
