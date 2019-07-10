using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Length
{
    /// <summary>
    /// Operates on the variations of the stress intensity factors ΔKI and ΔKII, which should be obtained by applying the 
    /// stress variation Δσ (or displacement variation Δu) on the structure.
    /// </summary>
    public class ParisLawIncrement2D: ICrackGrowthLengthLaw2D
    {
        private readonly int cyclesIncrement;
        private readonly double parisLawConstant;
        private readonly double parisLawExponent;

        public ParisLawIncrement2D(double parisLawConstant, double parisLawExponent, int cyclesIncrement)
        {
            // TODO: Add checks for the parameters.
            this.parisLawConstant = parisLawConstant;
            this.parisLawExponent = parisLawExponent;
            this.cyclesIncrement = cyclesIncrement;
        }

        public double ComputeGrowthLength(double sif1, double sif2)
        {
            double dKeff = CalcEffectiveSifRegular(sif1, sif2);
            return parisLawConstant * FindElapsedCycles() * Math.Pow(dKeff, parisLawExponent);
        }

        private int FindElapsedCycles()
        {
            return cyclesIncrement;
            //throw new NotImplementedException(
            //    "There should be a reference to an object that tracks the current cycle, e.g. the model");
        }

        private double CalcEffectiveSifRegular(double sif1, double sif2)
            => Math.Sqrt(sif1 * sif1 + sif2 * sif2);

        private double CalcEffectiveSifTanaka(double sif1, double sif2)
            => Math.Pow(Math.Pow(sif1, 4.0) + 8.0 * Math.Pow(sif2, 4.0), 0.25);
    }
}
