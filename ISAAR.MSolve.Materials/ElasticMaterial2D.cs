using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Materials
{
    public class ElasticMaterial2D : IIsotropicContinuumMaterial2D
    {
        private readonly double[] strains = new double[3];
        private readonly double[] stresses = new double[3];
        private double[,] constitutiveMatrix = null;

        public double[] Coordinates { get; set; }
        public double PoissonRatio { get; set; }
        public StressState2D StressState { get; }
        public double YoungModulus { get; set; }

        public ElasticMaterial2D(StressState2D stressState)
        {
            this.StressState = stressState;
        }

        #region IFiniteElementMaterial3D

        public ElasticityTensorContinuum2D ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new StressStrainVectorContinuum2D(new double[3]));
                return new ElasticityTensorContinuum2D(constitutiveMatrix);
            }
        }

        public StressStrainVectorContinuum2D Stresses { get { return new StressStrainVectorContinuum2D(stresses); } }

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

        public void UpdateMaterial(StressStrainVectorContinuum2D strains)
        {
            strains.CopyTo(this.strains, 0);
            constitutiveMatrix = new double[3, 3];
            if (StressState == StressState2D.PlaneStress)
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
        object ICloneable.Clone() => Clone();

        public ElasticMaterial2D Clone()
        {
            return new ElasticMaterial2D(StressState)
            {
                PoissonRatio = this.PoissonRatio,
                YoungModulus = this.YoungModulus
            };
        }

        #endregion

    }
}
