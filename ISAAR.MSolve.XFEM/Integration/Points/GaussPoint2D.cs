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
    class GaussPoint2D: INaturalPoint2D
    {
        public double Xi { get; }
        public double Eta { get; }
        public double Weight { get; }

        public GaussPoint2D(double xi, double eta, double weight)
        {
            this.Xi = xi;
            this.Eta = eta;
            this.Weight = weight;
        }
    }
}
