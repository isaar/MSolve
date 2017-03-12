using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    class NaturalPoint2D: INaturalPoint2D
    {
        public double Xi { get; }
        public double Eta { get; }

        public NaturalPoint2D(double xi, double eta)
        {
            this.Xi = xi;
            this.Eta = eta;
        }

        public override string ToString()
        {
            return "(" + Xi + " , " + Eta + ")";
        }
    }
}
