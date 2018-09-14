using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials.VonMisesStress
{
    /// <summary>
    /// Calculates the equivalent von Mises stress according to https://en.wikipedia.org/wiki/Von_Mises_yield_criterion#Summary.
    /// For plane strain conditions, the components s23 = s31 = 0, but s33 = E*v / ((1+v)*(1-2*v)) * (e11+e22). Therefore it
    /// depends on the Young modulus E and Poisson ratio v.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ElasticPlaneStrainVonMises: IVonMisesStress2D
    {
        private readonly double E;
        private readonly double v;

        public ElasticPlaneStrainVonMises(ElasticMaterial2D material)
        {
            this.E = material.YoungModulus;
            this.v = material.PoissonRatio;
        }

        public double Calculate(double[] strainTensor2D, double[] cauchyStressTensor2D)
        {
            double e11 = strainTensor2D[0];
            double e22 = strainTensor2D[1];
            double s11 = cauchyStressTensor2D[0];
            double s22 = cauchyStressTensor2D[1];
            double s12 = cauchyStressTensor2D[2];

            double s33 = E * v / ((1 + v) * (1 - 2 * v)) * (e11 + e22);
            double s11_22 = s11 - s22;
            double s22_33 = s22 - s33;
            double s33_11 = s33 - s11;
            return Math.Sqrt(0.5 * (s11_22 * s11_22 + s22_33 * s22_33 + s33_11 * s33_11 + 6 * s12 * s12));
        }
    }
}
