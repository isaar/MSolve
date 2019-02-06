using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces; //using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.Materials
{
    public class ElasticMaterial3D_v2DefGrad : IContinuumMaterial3DDefGrad
    {
        private readonly double[] strains = new double[6];
        private readonly double[] stresses = new double[6];
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

        public Vector Stresses { get { return new Vector(stresses); } }
        
        public Matrix2D ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new double[9]);
                return new Matrix2D(constitutiveMatrix);
            }
        }

        public void UpdateMaterial(double[] DefGradVec )
        {
            //throw new NotImplementedException();

            double[,] DGtr = new double[3, 3] { { DefGradVec [0], DefGradVec[7], DefGradVec[5] },
                                                { DefGradVec [3], DefGradVec[1], DefGradVec[8] },
                                                { DefGradVec [6], DefGradVec[4], DefGradVec[2] }};

            double[,] GL = new double[3, 3];
            double[] GLvec = new double[6];
            for(int m = 0; m < 3; m++)
            { 
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0;
                    for (int p = 0; p < 3; p++)
                    {
                        GL[m, n] += DGtr[m, p] * DGtr[n, p];
                    }
                }
            }
            for (int m = 0; m < 3; m++)
            {
                GL[m, m] += -1;
            }
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0.5 * GL[m, n];
                }
            }

            //
            for (int m = 0; m < 3; m++)
            {
                GLvec[m] = GL[m, m];
            }
            GLvec[3] = 2 * GL[0, 1];
            GLvec[4] = 2 * GL[1, 2];
            GLvec[5] = 2 * GL[2, 0];


            double[] strains = GLvec;
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
            return new ElasticMaterial3D_v2DefGrad() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }

        object ICloneable.Clone()
        {
            return Clone();
        }

        #endregion

    }
}
