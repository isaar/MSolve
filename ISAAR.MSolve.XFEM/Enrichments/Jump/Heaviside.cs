using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Enrichments.Jump
{
    class Heaviside : IJumpEnrichment
    {
        public double ValueAt(double signedDistance)
        {
            return (signedDistance >= 0) ? 1.0 : 0.0; // What should H(0) be? For now H(0) = 1.
        }

        public double DerivativeAt(double signedDistance)
        {
            return (signedDistance == 0) ? Double.PositiveInfinity : 0.0;
        }
    }
}
