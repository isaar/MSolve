using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Integration.Points
{
    /// <summary>
    /// Must be immutable to use across many elements that employ the same integration rule
    /// </summary>
    class GaussPoint : NaturalPoint
    {
        public double Weight { get; }

        public GaussPoint(double xi, double weight) : base(xi)
        {
            this.Weight = weight;
        }

        public GaussPoint(double xi, double eta, double weight): base(xi, eta)
        {
            this.Weight = weight;
        }
    }
}
