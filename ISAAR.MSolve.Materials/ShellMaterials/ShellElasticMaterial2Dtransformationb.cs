using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Materials
{
    /// <summary>
    /// Isotropic elastic shell material that accounts for no out of plane shear deformation
    /// Authors Gerasimos Sotiropoulos
    /// </summary>
    public class ShellElasticMaterial2Dtransformationb : IShellMaterial
    {
		public double[] NormalVectorV3 { get; set; }
		public double[] TangentVectorV1 { get; set; }
		public double[] TangentVectorV2 { get; set; }
		public double YoungModulus { get; set; }
		public double PoissonRatio { get; set; }
        Matrix transformationMatrix; // gia to shell

        private bool modified; 
		private double[,] CartesianConstitutiveMatrix;
		private double[] CartesianStresses = new double[6];

		object ICloneable.Clone() => Clone();

		public IShellMaterial Clone()
		{
			return new ShellElasticMaterial2Dtransformationb()
			{
				YoungModulus = this.YoungModulus,
				PoissonRatio = this.PoissonRatio,
			};
		}

		public void UpdateMaterial(double[] cartesianStrains) //TODO: rename cartesian strains to strains 
		{
			if (CartesianConstitutiveMatrix == null)
			{
				this.CalculateConstitutiveMatrix(Vector.CreateFromArray(TangentVectorV1), Vector.CreateFromArray(TangentVectorV2));
			}

			for (int l = 0; l < 3; l++)
			{
				CartesianStresses[l] = 0;
				for (int m = 0; m < 3; m++)
				{
					CartesianStresses[l] += CartesianConstitutiveMatrix[l, m] * cartesianStrains[m];
				}
			}
		}

		private void CalculateConstitutiveMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
		{
            this.CalculateTransformationMatrix(Vector.CreateFromArray(TangentVectorV1), Vector.CreateFromArray(TangentVectorV2));

            var OriginalConstitutiveMatrix = new double[3, 3];
            //if (StressState == StressState2D.PlaneStress)
            {
                double aux = YoungModulus / (1 - PoissonRatio * PoissonRatio);
                OriginalConstitutiveMatrix[0, 0] = aux;
                OriginalConstitutiveMatrix[1, 1] = aux;
                OriginalConstitutiveMatrix[0, 1] = PoissonRatio * aux;
                OriginalConstitutiveMatrix[1, 0] = PoissonRatio * aux;
                OriginalConstitutiveMatrix[2, 2] = (1 - PoissonRatio) / 2 * aux;
            }



            //         var auxMatrix1 = new Matrix2D(2, 2);
            //auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            //auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            //auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            //auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            //(Matrix2D inverse, double det) = auxMatrix1.Invert2x2AndDeterminant();

            //var constitutiveMatrix = new Matrix2D(new double[3, 3]
            //{
            //	{
            //		inverse[0,0]*inverse[0,0],
            //		this.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-this.PoissonRatio)*inverse[1,0]*inverse[1,0],
            //		inverse[0,0]*inverse[1,0]
            //	},
            //	{
            //		this.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-this.PoissonRatio)*inverse[1,0]*inverse[1,0],
            //		inverse[1,1]*inverse[1,1],
            //		inverse[1,1]*inverse[1,0]
            //	},
            //	{
            //		inverse[0,0]*inverse[1,0],
            //		inverse[1,1]*inverse[1,0],
            //		0.5*(1-this.PoissonRatio)*inverse[0,0]*inverse[1,1]+(1+this.PoissonRatio)*inverse[1,0]*inverse[1,0]
            //	},
            //});

            //         // Integrate over thickness takes into account multiplication *t but not (E/(1-(ni^2)) and it will be added here
            //         ConstitutiveMatrix.Scale(YoungModulus / (1 - Math.Pow(PoissonRatio, 2)));
            var constitutiveMatrix = (transformationMatrix.Transpose() * (Matrix.CreateFromArray(OriginalConstitutiveMatrix)) * transformationMatrix);
            CartesianConstitutiveMatrix = constitutiveMatrix.CopyToArray2D();
		}

        private void CalculateTransformationMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
        {
            var auxMatrix1 = Matrix.CreateZero(2, 2);  //auxMatrix: covariant metric coefficients gab
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            Matrix inverse = auxMatrix1.Invert(); //inverse: contravariant metric coefficients g_ab (ekthetis ta a,b)
                                                  //TODO: auxMatrix1.Invert2x2AndDeterminant(1e-20) for bad geometry

            //Contravariant base vectors
            double[][] G_i = new double[2][];
            for (int i1 = 0; i1 < 2; i1++)
            {
                G_i[i1] = new double[3];
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i1][i2] = inverse[i1, 0] * surfaceBasisVector1[i2] + inverse[i1, 1] * surfaceBasisVector2[i2];
                }
            }

            //Normalised covariant base vectors
            double[][] Ei = new double[2][];// to trito den xreiazetai

            Ei[0] =surfaceBasisVector1.CopyToArray();
            double G1_norm = surfaceBasisVector1.Norm2();
            for (int i1 = 0; i1 < 3; i1++) { Ei[0][i1] = Ei[0][i1] / G1_norm; }

            double G2_dot_E1 = 0;
            for (int i1 = 0; i1 < 3; i1++) { G2_dot_E1 += surfaceBasisVector2[i1] * Ei[0][i1]; }

            double[] projection = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { projection[i1] = G2_dot_E1 * Ei[0][i1]; }

            Ei[1] = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = surfaceBasisVector2[i1] - projection[i1]; }
            double norm1 = (Vector.CreateFromArray(Ei[1])).Norm2();
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = Ei[1][i1] / norm1; }

            double[,] EiDOTG_j = new double[2, 2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    EiDOTG_j[i1, i2] = Vector.CreateFromArray(Ei[i1]).DotProduct(Vector.CreateFromArray(G_i[i2]));
                }
            }

            transformationMatrix = Matrix.CreateFromArray(new double[3, 3] { {EiDOTG_j[0,0]*EiDOTG_j[0,0],EiDOTG_j[0,1]*EiDOTG_j[0,1],EiDOTG_j[0,0]*EiDOTG_j[0,1]  },
                 {EiDOTG_j[1,0]*EiDOTG_j[1,0],EiDOTG_j[1,1]*EiDOTG_j[1,1],EiDOTG_j[1,0]*EiDOTG_j[1,1]  },
                {2*EiDOTG_j[1,0]*EiDOTG_j[0,0],2*EiDOTG_j[1,1]*EiDOTG_j[0,1],EiDOTG_j[1,0]*EiDOTG_j[0,1]+EiDOTG_j[1,1]*EiDOTG_j[0,0]   } });
        }

        private bool CheckIfConstitutiveMatrixChanged()
		{
			return false;
		}

		public double[] Stresses 
		{
			get { return CartesianStresses; }
		}

		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (CartesianConstitutiveMatrix == null) UpdateMaterial(new double[6]);
				return Matrix.CreateFromArray(CartesianConstitutiveMatrix);
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

		public void ClearState()
		{
		}
		public void ClearStresses()
		{

		}

		public double[] Coordinates
		{

			get { throw new NotImplementedException(); }
			set { throw new InvalidOperationException(); }
		}

	}
}
