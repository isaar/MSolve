using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class VTKPoint: ICartesianPoint2D
    {
        public double X { get; }
        public double Y { get; }
        public double[] Coordinates { get { return new double[] { X, Y }; } }
        public int ID { get; }

        public VTKPoint(int id, ICartesianPoint2D originalNode)
        {
            this.ID = id;
            this.X = originalNode.X;
            this.Y = originalNode.Y;
        }
    }
}
