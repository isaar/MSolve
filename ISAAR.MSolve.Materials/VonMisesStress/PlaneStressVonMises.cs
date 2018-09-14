using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials.VonMisesStress
{
    /// <summary>
    /// Calculates the equivalent von Mises stress according to https://en.wikipedia.org/wiki/Von_Mises_yield_criterion#Summary.
    /// For plane stress conditions, the computation doesn't depend on the material properties.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PlaneStressVonMises : IVonMisesStress2D
    {
        public double Calculate(double[] strainTensor2D, double[] cauchyStressTensor2D)
        {
            double s11 = cauchyStressTensor2D[0];
            double s22 = cauchyStressTensor2D[1];
            double s12 = cauchyStressTensor2D[2];

            return Math.Sqrt(s11 * s11 - s11 * s22 + s22 * s22 + 3 * s12 * s12);
        }
    }
}
