using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Integration.Points
{
    class GaussPoint1D: INaturalPoint1D
    {
        public double Xi { get; }
        public double Weight { get; }

        public GaussPoint1D(double xi, double weight)
        {
            this.Xi = xi;
            this.Weight = weight;
        }
    }
}
