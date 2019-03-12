using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.MaterialInterpolation
{
    public class PowerLawMaterialInterpolation : IMaterialInterpolation
    {
        private readonly double penaltyExponent;
        private readonly double youngModulus;

        public PowerLawMaterialInterpolation(double youngModulus, double penaltyExponent, double minDensity)
        {
            this.youngModulus = youngModulus;
            this.penaltyExponent = penaltyExponent;
            this.MinDensity = minDensity;
        }

        public double MaxDensity => 1.0;

        public double MinDensity { get; }

        public double CalcMaterialProperty(double elementDensity)
            => Math.Pow(elementDensity, penaltyExponent) * youngModulus;

        public double CalcMaterialPropertyDerivative(double elementDensity)
            => penaltyExponent * Math.Pow(elementDensity, penaltyExponent - 1.0) * youngModulus;
    }
}
