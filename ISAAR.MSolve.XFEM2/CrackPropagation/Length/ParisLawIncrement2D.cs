using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Length
{
    class ParisLawIncrement2D: ICrackGrowthLengthLaw2D
    {
        private readonly double parisLawConstant;
        private readonly double parisLawExponent;

        public ParisLawIncrement2D(double parisLawConstant, double parisLawExponent)
        {
            // TODO: Add checks for the parameters.
            this.parisLawConstant = parisLawConstant;
            this.parisLawExponent = parisLawExponent;
        }

        public double ComputeGrowthLength(double sif1, double sif2)
        {
            double mixedModeSifCorrectionTanaka = Math.Pow(Math.Pow(sif1, 4.0) + 8.0 * Math.Pow(sif2, 4.0), 0.25);
            return parisLawConstant * FindElapsedCycles() * Math.Pow(mixedModeSifCorrectionTanaka, parisLawExponent);
        }

        private int FindElapsedCycles()
        {
            throw new NotImplementedException(
                "There should be a reference to an object that tracks the current cycle, e.g. the model");
        }
    }
}
