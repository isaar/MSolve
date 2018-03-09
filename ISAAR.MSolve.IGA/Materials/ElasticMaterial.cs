

using System;
using ISAAR.MSolve.IGA.Interfaces;

namespace ISAAR.MSolve.IGA.Problems.Structural.Constitutive
{
    public class ElasticMaterial : IIsogeometricMaterial
    {
        public double YoungModulus { get; set; }
        public double PoissonRatio { get; set; }
        public double[] Coordinates { get; set; }

        #region IFiniteElementMaterial Members

        public int ID
        {
            get { return 1; }
        }

        public bool Modified
        {
            get { return false; }
        }

        #endregion

        #region ICloneable Members

        public object Clone()
        {
            return new ElasticMaterial() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }

        public void ResetModified()
        {
        }

        #endregion
    }
}
