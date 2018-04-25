using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Materials
{
    public class ElasticMaterial2D : IFiniteElementMaterial3D
    {
        private readonly double[] strains = new double[3];
        private readonly double[] stresses = new double[3];
        private double [,] constitutiveMatrix = null;
        public double YoungModulus { get; set; }
        public double PoissonRatio { get; set; }
        public String StressState { get; set; }
        public double[] Coordinates { get; set; }


        #region IFiniteElementMaterial3D

        public IMatrix2D ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new double[3]);
                return new Matrix2D(constitutiveMatrix);
            }
        }

        public double[] Stresses { get { return stresses; } }

        public void ClearState()
        {
            throw new NotImplementedException();
        }

        public void ClearStresses()
        {
            throw new NotImplementedException();
        }

        public void SaveState()
        {
            throw new NotImplementedException();
        }

        public void UpdateMaterial(double[] strains)
        {
            strains.CopyTo(this.strains, 0);
            constitutiveMatrix = new double[3, 3];
            if (StressState.ToLower().Equals("plstress"))
            {
                double aux = YoungModulus / (1 - PoissonRatio * PoissonRatio);
                constitutiveMatrix[0, 0] = aux;
                constitutiveMatrix[1, 1] = aux;
                constitutiveMatrix[0, 1] = PoissonRatio * aux;
                constitutiveMatrix[1, 0] = PoissonRatio * aux;
                constitutiveMatrix[2, 2] = (1 - PoissonRatio) / 2 * aux;
            }
            else
            {
                double aux = YoungModulus / (1 + PoissonRatio) / (1 - 2 * PoissonRatio);
                constitutiveMatrix[0, 0] = aux * (1 - PoissonRatio);
                constitutiveMatrix[1, 1] = aux * (1 - PoissonRatio);
                constitutiveMatrix[0, 1] = PoissonRatio * aux;
                constitutiveMatrix[1, 0] = PoissonRatio * aux;
                constitutiveMatrix[2, 2] = (1 - 2 * PoissonRatio) / 2 * aux;
            }
        }

        #endregion

        #region IFiniteElementMaterial

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
            return new ElasticMaterial2D() {YoungModulus=this.YoungModulus, PoissonRatio =this.PoissonRatio, StressState=this.StressState };
        }

        #endregion

    }
}
