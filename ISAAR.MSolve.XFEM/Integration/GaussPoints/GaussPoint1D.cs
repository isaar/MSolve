using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Integration.GaussPoints
{
    class GaussPoint1D
    {
        public double X { get; }
        public double Weight { get; }

        public GaussPoint1D(double x, double weight)
        {
            this.X = x;
            this.Weight = weight;
        }
    }
}
