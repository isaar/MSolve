using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Enrichments.Jump
{
    class PolynomialSmoothedHeaviside : IJumpEnrichment
    {
        private readonly double e;

        public PolynomialSmoothedHeaviside(double e)
        {
            if (e <= 0.0)
            {
                throw new ArgumentException("The half length of the smoothed region must be e > 0 but was e = " + e);
            }
            this.e = e;
        }

        public double ValueAt(double signedDistance)
        {
            if (signedDistance <= -e) return 0.0;
            else if (signedDistance >= e) return 1.0;
            else return 0.5 + 0.75 * 
                    (signedDistance / e - (signedDistance * signedDistance * signedDistance) / (3.0 * e * e * e ));
        }

        public double DerivativeAt(double signedDistance)
        {
            if (signedDistance <= -e || signedDistance >= e) return 0.0;
            else return 0.75 * (e - signedDistance * signedDistance / e);
        }
    }
}
