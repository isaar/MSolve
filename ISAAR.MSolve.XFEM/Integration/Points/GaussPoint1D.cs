using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Integration.Points
{
    /// <summary>
    /// Must be immutable to use across many elements that employ the same integration rule
    /// </summary>
    class GaussPoint1D : INaturalPoint1D
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
