using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;

namespace ISAAR.MSolve.Materials
{
    /// <summary>
    /// Isotropic.
    /// Authors Gerasimos-Serafeim 
    /// </summary>
    public class ShellElasticMaterial3D : IShellMaterial 
    {
        public double[] NormalVectorV3 { get; set; }
        public double[] TangentVectorV1 { get; set; }
	    public double[] TangentVectorV2 { get; set; }
		public double YoungModulus { get; set; }
        public double PoissonRatio { get; set; }
        public double ShearCorrectionCoefficientK { get; set; }

        private bool modified; // opws sto MohrCoulomb gia to modified
        private double [,] ConsCartes;
        private double[] SPKvec = new double[6];

        object ICloneable.Clone() => Clone();

        public IShellMaterial Clone()
        {
            return new ShellElasticMaterial3D()
            {
               YoungModulus=this.YoungModulus,
               PoissonRatio=this.PoissonRatio,
               ShearCorrectionCoefficientK= this.ShearCorrectionCoefficientK              
            };
        }

        public void UpdateMaterial(double[] GLvec)
        {
            //elegxos initialize matrices
            if (ConsCartes==null)
            {
                this.CalculateConstituveMatrixInCartesianCoordinates(NormalVectorV3,TangentVectorV1);                
            }

            for (int l = 0; l < 6; l++)
            {
                SPKvec[l] = 0;
                for (int m = 0; m < 6; m++)
                {
                    SPKvec[l] += ConsCartes[l, m] * GLvec[m];
                }
            }
        }

        private bool CheckIfConstitutiveMatrixChanged()
        {
            return false;
        }

        public double[] Stresses // opws xrhsimopoeitai sto mohrcoulomb kai hexa8
        {
            get { return SPKvec; }
        }

        public IMatrix2D  ConstitutiveMatrix
        {
            get
            {
                if (ConsCartes == null) UpdateMaterial(new double[6]);
                return new Matrix2D(ConsCartes);
            } 
        }

        public void SaveState()
        {           
        }

        public bool Modified => modified;

        public void ResetModified()
        {
            modified = false;
        }

        public int ID
        {
            get { throw new NotImplementedException(); }
        }

        public void ClearState() // pithanws TODO 
        {
        }
        public void ClearStresses()
        {

        }

        public double[] Coordinates
        {

            get { throw new NotImplementedException();  }
            set { throw new InvalidOperationException(); }
        }

        private void CalculateConstituveMatrixInCartesianCoordinates(double[] V3, double[] V1)
        {
            double E = YoungModulus;
            double ni = PoissonRatio;
            ConsCartes = new double[6, 6];
            double[,] Cons = new double[6, 6];
            double[,] Cons_T_e = new double[6, 6];
            double[] V2 = new double[3];



            for (int k = 0; k < 2; k++) { Cons[k, k] = E / (1 - Math.Pow(ni, 2)); }
            Cons[0, 1] = ni * E / (1 - Math.Pow(ni, 2));
            Cons[1, 0] = ni * E / (1 - Math.Pow(ni, 2));
            Cons[3, 3] = (1 - ni) * (0.5) * E / (1 - Math.Pow(ni, 2));
            Cons[4, 4] = (1 - ni) * (0.5) * E / (1 - Math.Pow(ni, 2)); //Cons[4, 4] = (1 - ni) * (0.41666666667) * E / (1 - Math.Pow(ni, 2));
            Cons[5, 5] = (1 - ni) * (0.5) * E / (1 - Math.Pow(ni, 2)); //Cons[5, 5] = (1 - ni) * (0.41666666667) * E / (1 - Math.Pow(ni, 2));
            //for (int k = 0; k < 2; k++)
            //{ Cons[4 + k, 4 + k] = (5 / 6) * (1 - ni) * (0.5) * E / (1 - Math.Pow(ni, 2)); }

            V2[0] = V3[1] * V1[2] - V3[2] * V1[1];
            V2[1] = V3[2] * V1[0] - V3[0] * V1[2];
            V2[2] = V3[0] * V1[1] - V3[1] * V1[0];

            double l1 = V1[0];
            double m1 = V1[1];
            double n1 = V1[2];

            double l2 = V2[0];
            double m2 = V2[1];
            double n2 = V2[2];

            double l3 = V3[0];
            double m3 = V3[1];
            double n3 = V3[2];

            double[,] T_e = new double[6, 6];
            for (int i = 0; i < 3; i++)
            {
                T_e[0, i] = (V1[i] * V1[i]);
                T_e[1, i] = (V2[i] * V2[i]);
                T_e[2, i] = (V3[i] * V3[i]);

                T_e[3, i] = (2 * V1[i] * V2[i]);
                T_e[4, i] = (2 * V2[i] * V3[i]);
                T_e[5, i] = (2 * V3[i] * V1[i]);

                T_e[0, 3 + i] = (V1[i] * V1[1 + i - 3 * i * (i - 1) / 2]);
                T_e[1, 3 + i] = (V2[i] * V2[1 + i - 3 * i * (i - 1) / 2]);
                T_e[2, 3 + i] = (V3[i] * V3[1 + i - 3 * i * (i - 1) / 2]);

                T_e[3, 3 + i] = (V1[i] * V2[1 + i - 3 * i * (i - 1) / 2] + V2[i] * V1[1 + i - 3 * i * (i - 1) / 2]);
                T_e[4, 3 + i] = (V2[i] * V3[1 + i - 3 * i * (i - 1) / 2] + V3[i] * V2[1 + i - 3 * i * (i - 1) / 2]);
                T_e[5, 3 + i] = (V3[i] * V1[1 + i - 3 * i * (i - 1) / 2] + V1[i] * V3[1 + i - 3 * i * (i - 1) / 2]);
            }

            // multiplication [Te']*[cons]*[Te];

            for (int i = 0; i < 6; i++) //TODO use if LinearAlgebra
            {
                for (int k = 0; k < 6; k++)
                {
                    Cons_T_e[i, k] = 0;
                    for (int l = 0; l < 6; l++)
                    { Cons_T_e[i, k] += Cons[i, l] * T_e[l, k]; }
                }
            }

            for (int i = 0; i < 6; i++)
            {
                for (int k = 0; k < 6; k++)
                {
                    ConsCartes[i, k] = 0;
                    for (int l = 0; l < 6; l++)
                    { ConsCartes[i, k] += T_e[l, i] * Cons_T_e[l, k]; }
                }
            }
        }
    }
}
