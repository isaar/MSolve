using System;
using System.Collections.Generic;
using System.Text;

//TODO: it will probably need to be removed if it turns out that the material properties at each integration point are needed for 
//      plane strain. Then the von Mises stress should be queried from the material of each integration point.
namespace ISAAR.MSolve.Materials.VonMisesStress
{
    /// <summary>
    /// Calculates the von Mises equivalent stress given the Cauchy stress.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IVonMisesStress2D
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="strainTensor2D">Its components are: e11, e22, e12</param>
        /// <param name="cauchyStressTensor2D">Its components are: s11, s22, s12.</param>
        /// <returns></returns>
        double Calculate(double[] strainTensor2D, double[] cauchyStressTensor2D);
    }
}
