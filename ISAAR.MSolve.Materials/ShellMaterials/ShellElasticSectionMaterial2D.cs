using System;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.Materials
{
    public class ShellElasticSectionMaterial2D : IShellSectionMaterial
	{
		object ICloneable.Clone()
		{
			return Clone();
		}

		public double[] MembraneForces { get; }
		public double[] Moments { get; }
		public IMatrixView MembraneConstitutiveMatrix { get; private set; }
		public IMatrixView BendingConstitutiveMatrix { get; private set; }
		public IMatrixView CouplingConstitutiveMatrix { get; private set; }
		


		public double[] NormalVectorV3 { get; set; }
		public double[] TangentVectorV1 { get; set; }
		public double[] TangentVectorV2 { get; set; }
		public double Thickness { get; set; }

		public void UpdateMaterial(double[] membraneStrains, double[] bendingStrains)
		{
			if (MembraneConstitutiveMatrix == null)
			{
				CalculateMembraneConstitutiveMatrix(TangentVectorV1, TangentVectorV2, Thickness);
			}

			if (BendingConstitutiveMatrix == null)
			{
				CalculateBendingConstitutiveMatrix(TangentVectorV1, TangentVectorV2, Thickness);
			}

			if (CouplingConstitutiveMatrix == null)
			{
				CalculateCouplingConstitutiveMatrix(TangentVectorV1, TangentVectorV2, Thickness);
			}

			for (int l = 0; l < 3; l++)
			{
				MembraneForces[l] = 0;
				for (int m = 0; m < 3; m++)
				{
					MembraneForces[l] += MembraneConstitutiveMatrix[l, m] * membraneStrains[m];
				}
			}

			for (int l = 0; l < 3; l++)
			{
				Moments[l] = 0;
				for (int m = 0; m < 3; m++)
				{
					Moments[l] += BendingConstitutiveMatrix[l, m] * bendingStrains[m];
				}
			}
		}

		private void CalculateCouplingConstitutiveMatrix(double[] vector1, double[] vector2, double thickness)
		{
			CouplingConstitutiveMatrix = Matrix.CreateZero(3,3);
		}

		private void CalculateBendingConstitutiveMatrix(double[] vector1, double[] vector2, double thickness)
		{
			BendingConstitutiveMatrix = CalculateConstitutiveMatrix(vector1, vector2);
			BendingConstitutiveMatrix.Scale(YoungModulus * Math.Pow(Thickness, 3) /
			                                12 / (1 - Math.Pow(PoissonRatio, 2)));
		}

		private void CalculateMembraneConstitutiveMatrix(double[] vector1, double[] vector2, double thickness)
		{
			MembraneConstitutiveMatrix= CalculateConstitutiveMatrix(vector1,vector2);
			MembraneConstitutiveMatrix.Scale(YoungModulus * Thickness /
			                                 (1 - Math.Pow(PoissonRatio, 2)));
		}

		private Matrix CalculateConstitutiveMatrix(double[] surfaceBasisVector1, double[] surfaceBasisVector2)
		{
			var auxMatrix1 = Matrix.CreateZero(2, 2);
			auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
			auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
			auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
			auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            (Matrix inverse, double det) = auxMatrix1.InvertAndDeterminant();

			var constitutiveMatrix = Matrix.CreateFromArray(new double[3, 3]
			{
				{
					inverse[0,0]*inverse[0,0],
					this.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-this.PoissonRatio)*inverse[1,0]*inverse[1,0],
					inverse[0,0]*inverse[1,0]
				},
				{
					this.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-this.PoissonRatio)*inverse[1,0]*inverse[1,0],
					inverse[1,1]*inverse[1,1],
					inverse[1,1]*inverse[1,0]
				},
				{
					inverse[0,0]*inverse[1,0],
					inverse[1,1]*inverse[1,0],
					0.5*(1-this.PoissonRatio)*inverse[0,0]*inverse[1,1]+(1+this.PoissonRatio)*inverse[1,0]*inverse[1,0]
				},
			});
			return constitutiveMatrix;
		}


		public int ID { get; }
		public bool Modified { get; private set; }

		public void ResetModified()
		{
			Modified = false;
		}

		public double[] Coordinates { get; set; }
		public double YoungModulus { get; set; }
		public double PoissonRatio { get; set; }
		public void SaveState()
		{
		}

		public void ClearState()
		{
		}

		public void ClearStresses()
		{
		}

		public IShellSectionMaterial Clone()
		{
			return new ShellElasticSectionMaterial2D()
			{
				YoungModulus = this.YoungModulus,
				PoissonRatio = this.PoissonRatio,
				Thickness = this.Thickness
			};

		}
	}
}
