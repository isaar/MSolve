using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.Materials
{
	 public class ShellElasticMaterial2D : IShellMaterial
	{
		public double[] NormalVectorV3 { get; set; }
		public double[] TangentVectorV1 { get; set; }
		public double[] TangentVectorV2 { get; set; }
		public double YoungModulus { get; set; }
		public double PoissonRatio { get; set; }

		private bool modified; 
		private Matrix CartesianConstitutiveMatrix;
		private double[] CartesianStresses = new double[6];

		object ICloneable.Clone() => Clone();

		public IShellMaterial Clone()
		{
			return new ShellElasticMaterial2D()
			{
				YoungModulus = this.YoungModulus,
				PoissonRatio = this.PoissonRatio,
			};
		}

		public void UpdateMaterial(double[] cartesianStrains)
		{
			if (CartesianConstitutiveMatrix == null)
			{
				this.CalculateConstitutiveMatrix(TangentVectorV1, TangentVectorV2);
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

		private void CalculateConstitutiveMatrix(double[] surfaceBasisVector1, double[] surfaceBasisVector2)
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
					PoissonRatio*inverse[0,0]*inverse[1,1]+(1-PoissonRatio)*inverse[1,0]*inverse[1,0],
					inverse[0,0]*inverse[1,0]
				},
				{
					PoissonRatio*inverse[0,0]*inverse[1,1]+(1-PoissonRatio)*inverse[1,0]*inverse[1,0],
					inverse[1,1]*inverse[1,1],
					inverse[1,1]*inverse[1,0]
				},
				{
					inverse[0,0]*inverse[1,0],
					inverse[1,1]*inverse[1,0],
					0.5*(1-PoissonRatio)*inverse[0,0]*inverse[1,1]+(1+PoissonRatio)*inverse[1,0]*inverse[1,0]
				},
			});
			constitutiveMatrix.ScaleIntoThis(YoungModulus/(1-Math.Pow(PoissonRatio,2)));
			CartesianConstitutiveMatrix = constitutiveMatrix;
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
				return CartesianConstitutiveMatrix;
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
