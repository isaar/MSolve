using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    class PolarPoint2D
    {
        public double R { get; }
        public double Theta { get; }

        public PolarPoint2D (double r, double theta)
        {
            this.R = r;
            this.Theta = theta;
        }
    }
}
