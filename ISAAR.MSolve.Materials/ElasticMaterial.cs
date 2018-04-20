using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Materials
{
    public class ElasticMaterial : IFiniteElementMaterial
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

        public void ResetModified()
        {
        }

        #endregion

        #region ICloneable Members

        public object Clone()
        {
            return new ElasticMaterial() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }

        #endregion
    }
}
