using System;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.Materials
{
    public class ElasticMaterial3D : IIsotropicContinuumMaterial3D
    {
        //private readonly double[] strains = new double[6];
        private readonly double[] stresses = new double[6];
        private Matrix constitutiveMatrix = null;
        public double YoungModulus { get; set; }
        public double PoissonRatio { get; set; }
        public double[] Coordinates { get; set; }
        private readonly double[] incrementalStrains = new double[6];
        private double[] stressesNew = new double[6];

        private Matrix GetConstitutiveMatrix()
        {
            double fE1 = YoungModulus / (double)(1 + PoissonRatio);
            double fE2 = fE1 * PoissonRatio / (double)(1 - 2 * PoissonRatio);
            double fE3 = fE1 + fE2;
            double fE4 = fE1 * 0.5;
            var afE = Matrix.CreateZero(6, 6);
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

            return afE;
        }

        private void CalculateNextStressStrainPoint()
        {
            var stressesElastic = new double[6];
            for (int i = 0; i < 6; i++)
            {
                stressesElastic[i] = this.stresses[i];
                for (int j = 0; j < 6; j++)
                    stressesElastic[i] += this.constitutiveMatrix[i, j] * this.incrementalStrains[j];
            }

            this.stressesNew = stressesElastic;
        }

        #region IFiniteElementMaterial Members

        public int ID => 1;

        public bool Modified => false;

        public void ResetModified() { }

        #endregion

        #region IFiniteElementMaterial3D Members

        public double[] Stresses => stressesNew;

        public IMatrixView ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new double[6]);
                return constitutiveMatrix;
            }
        }

        public void UpdateMaterial(double[] strainsIncrement)
        {
            //throw new NotImplementedException();
            this.incrementalStrains.CopyFrom(strainsIncrement);
            constitutiveMatrix = GetConstitutiveMatrix();
            this.CalculateNextStressStrainPoint();

        }

        public void ClearState()
        {
            //constitutiveMatrix.Clear();
            incrementalStrains.Clear();
            stresses.Clear();
            stressesNew.Clear();
        }

        public void SaveState() => stresses.CopyFrom(stressesNew);

        public void ClearStresses()
        {
            stresses.Clear();
            stressesNew.Clear();
        }

        #endregion

        #region ICloneable Members

        object ICloneable.Clone() => Clone();

        public ElasticMaterial3D Clone()
        {
            return new ElasticMaterial3D() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }

        #endregion

    }

}
