using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Materials
{
    public class ElasticMaterial : IFiniteElementMaterial
    {
        private readonly double[] strains = new double[3];
        private readonly double[] incrementalStrains = new double[3];
        private readonly double[] stresses = new double[3];
        private double[] stressesNew = new double[3];
        private double[,] constitutiveMatrix = null;
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

        public ElasticMaterial Clone()
        {
            return new ElasticMaterial() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }

        public void SaveState()
        {
            Array.Copy(this.stressesNew, this.stresses, 3);
        }

        public void ClearState()
        {
            if (constitutiveMatrix != null) Array.Clear(constitutiveMatrix, 0, constitutiveMatrix.Length);
            Array.Clear(incrementalStrains, 0, incrementalStrains.Length);
            Array.Clear(stresses, 0, stresses.Length);
            Array.Clear(stressesNew, 0, stressesNew.Length);
        }

        public void ClearStresses()
        {
            throw new NotImplementedException();
        }

        #endregion
    }
}
