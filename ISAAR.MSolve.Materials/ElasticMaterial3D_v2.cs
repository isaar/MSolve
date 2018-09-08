using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; 
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces; 
using ISAAR.MSolve.Numerical.LinearAlgebra; 

namespace ISAAR.MSolve.Materials
{
    public class ElasticMaterial3D_v2 : IContinuumMaterial3D, IIsotropicContinuumMaterial3D
    {
        private readonly double[] strains = new double[6];
        private readonly StressStrainVectorContinuum3D stresses = new StressStrainVectorContinuum3D();
        private double[,] constitutiveMatrix = null;
        public double YoungModulus { get; set; }
        public double PoissonRatio { get; set; }
        public double[] Coordinates { get; set; }

        private double[,] GetConstitutiveMatrix()
        {
            double fE1 = YoungModulus / (double)(1 + PoissonRatio);
            double fE2 = fE1 * PoissonRatio / (double)(1 - 2 * PoissonRatio);
            double fE3 = fE1 + fE2;
            double fE4 = fE1 * 0.5;
            double[,] afE = new double[6, 6];
            afE[0, 0] = fE3;
            afE[0, 1] = fE2;
            afE[0, 2] = fE2;
            afE[1, 0] = fE2;
            afE[1, 1] = fE3;
            afE[1, 2] = fE2;
            afE[2, 0] = fE2;
            afE[2, 1] = fE2;
            afE[2, 2] = fE3;
            afE[3, 3] = fE4;
            afE[4, 4] = fE4;
            afE[5, 5] = fE4;

            Vector s = (new Matrix2D(afE)) * (new Vector(strains));
            s.Data.CopyTo(stresses, 0);

            return afE;
        }

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

        #region IFiniteElementMaterial3D Members

        public StressStrainVectorContinuum3D Stresses { get { return stresses; } }

        public ElasticityTensorContinuum3D ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new double[6]);
                return new ElasticityTensorContinuum3D(constitutiveMatrix);
            }
        }

        public void UpdateMaterial(double[] strains)
        {
            //throw new NotImplementedException();

            strains.CopyTo(this.strains, 0);
            constitutiveMatrix = GetConstitutiveMatrix();
        }

        public void ClearState()
        {
            //throw new NotImplementedException();
        }

        public void SaveState()
        {
            //throw new NotImplementedException();
        }

        public void ClearStresses()
        {
            //throw new NotImplementedException();
        }

        #endregion

        #region ICloneable Members

        public object Clone()
        {
            return new ElasticMaterial3D_v2() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }

		public void UpdateMaterial(StressStrainVectorContinuum3D strains)
		{
			throw new NotImplementedException();
		}

		#endregion

	}
}
