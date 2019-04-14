using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.MaterialInterpolation
{
    public class ModifiedPowerLawMaterialInterpolation : IMaterialInterpolation
    {
        private readonly double penaltyExponent;
        private readonly double youngModulusSolid;
        private readonly double youngModulusVoid;

        public ModifiedPowerLawMaterialInterpolation(double youngModulusSolid, double youngModulusVoid, double penaltyExponent)
        {
            this.youngModulusSolid = youngModulusSolid;
            this.youngModulusVoid = youngModulusVoid;
            this.penaltyExponent = penaltyExponent;
        }

        public double MaxDensity => 1.0;

        public double MinDensity => 0.0;

        public double CalcMaterialProperty(double elementDensity)
            => youngModulusVoid + Math.Pow(elementDensity, penaltyExponent) * (youngModulusSolid - youngModulusVoid);

        public double CalcMaterialPropertyDerivative(double elementDensity)
            => penaltyExponent * Math.Pow(elementDensity, penaltyExponent - 1.0) * (youngModulusSolid - youngModulusVoid);
    }
}
