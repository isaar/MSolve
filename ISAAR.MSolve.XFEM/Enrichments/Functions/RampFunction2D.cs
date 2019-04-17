using System;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Functions
{
    class RampFunction2D: IEnrichmentFunction2D
    {
        public RampFunction2D()
        {
        }

        public double EvaluateAt(double signedDistance)
        {
            return Math.Abs(signedDistance);
        }

        public EvaluatedFunction2D EvaluateAllAt(double signedDistance, Vector2 normalVector)
        {
            int sign = Math.Sign(signedDistance);
            return new EvaluatedFunction2D(Math.Abs(signedDistance), sign * normalVector);
        }

        public override string ToString()
        {
            return "Ramp";
        }
    }
}
