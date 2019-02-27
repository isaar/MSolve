using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.IGA.Elements
{
	public class NURBSElement3D : Element, IStructuralIsogeometricElement
	{
		protected readonly static DOFType[] controlPointDOFTypes = new DOFType[] {DOFType.X, DOFType.Y, DOFType.Z};
		protected DOFType[][] dofTypes;
		protected IElementDofEnumerator_v2 dofEnumerator = new GenericDofEnumerator_v2();
		private DynamicMaterial dynamicProperties;


		public IElementDofEnumerator_v2 DofEnumerator
		{
			get { return dofEnumerator; }

			set { this.dofEnumerator = value; }
		}

		public ElementDimensions ElementDimensions
		{
			get { return ElementDimensions.ThreeD; }
		}

		public IList<IList<DOFType>> GetElementDOFTypes(IElement_v2 element)
		{
			var nurbsElement = (NURBSElement3D) element;
			dofTypes = new DOFType[nurbsElement.ControlPoints.Count][];
			for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
			{
				dofTypes[i] = controlPointDOFTypes;
			}

			return dofTypes;
		}

		public bool MaterialModified
		{
			get { throw new NotImplementedException(); }
		}

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements,
			double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		public void ClearMaterialState()
		{
			throw new NotImplementedException();
		}

		public IMatrix DampingMatrix(IElement_v2 element)
		{
			throw new NotImplementedException();
		}

		public IMatrix MassMatrix(IElement_v2 element)
		{
			throw new NotImplementedException();
		}

		//public IMatrix StiffnessMatrix(IElement_v2 element)
		//{
		//	var nurbsElement = (NURBSElement3D)element;
		//	var controlPoints = nurbsElement.ControlPoints.ToArray();
		//	IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(nurbsElement);
		//	Matrix2D stiffnessMatrixElement = new Matrix2D(controlPoints.Length * 3, controlPoints.Length * 3);

		//	NURBS3D nurbs = new NURBS3D(nurbsElement, controlPoints);

		//	for (int j = 0; j < gaussPoints.Count; j++)
		//	{
		//		var jacobianMatrix = CalculateJacobian(controlPoints, nurbs, j);

		//		double jacdet = CalculateJacobianDeterminant(jacobianMatrix);

		//		var inverseJacobian = CalculateInverseJacobian(jacobianMatrix, jacdet);

		//		var B1 = CalculateDeformationMatrix1(inverseJacobian);

		//		var B2 = CalculateDeformationMatrix2(controlPoints, nurbs, j);

		//		var B = new Matrix2D(B1) * new Matrix2D(B2);

		//		Matrix2D stiffnessMatrixGaussPoint = B.Transpose() * ((IContinuumMaterial3D)nurbsElement.Patch.Material).ConstitutiveMatrix * B * jacdet * gaussPoints[j].WeightFactor;

		//		for (int m = 0; m < controlPoints.Length * 3; m++)
		//		{
		//			for (int n = 0; n < controlPoints.Length * 3; n++)
		//			{
		//				stiffnessMatrixElement[m, n] += stiffnessMatrixGaussPoint[m, n];
		//			}
		//		}
		//	}
		//	return Matrix.CreateFromArray(stiffnessMatrixElement.Data);

		//}


		public IMatrix StiffnessMatrix(IElement_v2 element)
		{
			var nurbsElement = (NURBSElement3D)element;
			var controlPoints = nurbsElement.ControlPoints.ToArray();
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(nurbsElement);

			NURBS3D nurbs = new NURBS3D(nurbsElement, controlPoints);


			var bRows = 6;
			var bCols = controlPoints.Length * 3;
			var stiffnessMatrix = new double[bCols, bCols];
			var BTranspose = new double[bCols, bRows];
			var BTransposeMultStiffness = new double[bCols, bRows];

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(controlPoints, nurbs, j);

				double jacdet = CalculateJacobianDeterminant(jacobianMatrix);

				var inverseJacobian = CalculateInverseJacobian(jacobianMatrix, jacdet);

				var B1 = CalculateDeformationMatrix1(inverseJacobian);

				var B2 = CalculateDeformationMatrix2(controlPoints, nurbs, j);

				var B = new double[B1.GetLength(0), B2.GetLength(1)];

				for (int i = 0; i < B1.GetLength(0); i++)
					for (int k = 0; k < B2.GetLength(1); k++)
						for (int m = 0; m < B2.GetLength(0); m++)
							B[i, k] += B1[i, m] * B2[m, k];



				double wFactor = jacdet * gaussPoints[j].WeightFactor;

				Array.Clear(BTranspose, 0, bRows * bCols);

				for (int i = 0; i < bRows; i++)
					for (int k = 0; k < bCols; k++)
						BTranspose[k, i] = B[i, k] * wFactor;

				double tempcm = 0;
				double tempb = 0;
				var E = ((IContinuumMaterial3D)nurbsElement.Patch.Material).ConstitutiveMatrix.Data;
				Array.Clear(BTransposeMultStiffness, 0, bRows * bCols);
				for (int i = 0; i < bCols; i++)
				{
					for (int k = 0; k < bRows; k++)
					{
						tempb = BTranspose[i, k];
						for (int m = 0; m < bRows; m++)
						{
							tempcm = E[k, m];
							BTransposeMultStiffness[i, m] += tempb * tempcm;
						}
					}
				}

				tempb = 0;
				for (int i = 0; i < bCols; i++)
				{
					for (int k = 0; k < bRows; k++)
					{
						tempb = BTransposeMultStiffness[i, k];
						for (int m = 0; m < bCols; m++)
						{
							stiffnessMatrix[i, m] += tempb * B[k, m];
						}
					}
				}
			}

			return Matrix.CreateFromArray(stiffnessMatrix);
		}

		private static double[,] CalculateDeformationMatrix2(IList<ControlPoint> elementControlPoints, NURBS3D nurbs,
			int j)
		{
			double[,] B2 = new double[9, 3 * elementControlPoints.Count];
			for (int column = 0; column < 3 * elementControlPoints.Count; column += 3)
			{
				B2[0, column] += nurbs.NurbsDerivativeValuesKsi[column / 3, j];
				B2[1, column] += nurbs.NurbsDerivativeValuesHeta[column / 3, j];
				B2[2, column] += nurbs.NurbsDerivativeValuesZeta[column / 3, j];

				B2[3, column + 1] += nurbs.NurbsDerivativeValuesKsi[column / 3, j];
				B2[4, column + 1] += nurbs.NurbsDerivativeValuesHeta[column / 3, j];
				B2[5, column + 1] += nurbs.NurbsDerivativeValuesZeta[column / 3, j];

				B2[6, column + 2] += nurbs.NurbsDerivativeValuesKsi[column / 3, j];
				B2[7, column + 2] += nurbs.NurbsDerivativeValuesHeta[column / 3, j];
				B2[8, column + 2] += nurbs.NurbsDerivativeValuesZeta[column / 3, j];
			}

			return B2;
		}

		private static double[,] CalculateDeformationMatrix1(double[,] inverseJacobian)
		{
			var B1 = new double[6, 9];

			B1[0, 0] += inverseJacobian[0, 0];
			B1[0, 1] += inverseJacobian[0, 1];
			B1[0, 2] += inverseJacobian[0, 2];

			B1[1, 3] += inverseJacobian[1, 0];
			B1[1, 4] += inverseJacobian[1, 1];
			B1[1, 5] += inverseJacobian[1, 2];

			B1[2, 6] += inverseJacobian[2, 0];
			B1[2, 7] += inverseJacobian[2, 1];
			B1[2, 8] += inverseJacobian[2, 2];

			B1[3, 0] += inverseJacobian[1, 0];
			B1[3, 1] += inverseJacobian[1, 1];
			B1[3, 2] += inverseJacobian[1, 2];
			B1[3, 3] += inverseJacobian[0, 0];
			B1[3, 4] += inverseJacobian[0, 1];
			B1[3, 5] += inverseJacobian[0, 2];

			B1[4, 3] += inverseJacobian[2, 0];
			B1[4, 4] += inverseJacobian[2, 1];
			B1[4, 5] += inverseJacobian[2, 2];
			B1[4, 6] += inverseJacobian[1, 0];
			B1[4, 7] += inverseJacobian[1, 1];
			B1[4, 8] += inverseJacobian[1, 2];

			B1[5, 0] += inverseJacobian[2, 0];
			B1[5, 1] += inverseJacobian[2, 1];
			B1[5, 2] += inverseJacobian[2, 2];
			B1[5, 6] += inverseJacobian[0, 0];
			B1[5, 7] += inverseJacobian[0, 1];
			B1[5, 8] += inverseJacobian[0, 2];
			return B1;
		}

		private static double[,] CalculateInverseJacobian(double[,] jacobianMatrix, double jacdet)
		{
			var inverseJacobian = new double[3, 3];

			inverseJacobian[0, 0] =
				(jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[1, 2] * jacobianMatrix[2, 1]) / jacdet;
			inverseJacobian[0, 1] =
				(jacobianMatrix[0, 2] * jacobianMatrix[2, 1] - jacobianMatrix[0, 1] * jacobianMatrix[2, 2]) / jacdet;
			inverseJacobian[0, 2] =
				(jacobianMatrix[0, 1] * jacobianMatrix[1, 2] - jacobianMatrix[0, 2] * jacobianMatrix[1, 1]) / jacdet;

			inverseJacobian[1, 0] =
				(jacobianMatrix[1, 2] * jacobianMatrix[2, 0] - jacobianMatrix[1, 0] * jacobianMatrix[2, 2]) / jacdet;
			inverseJacobian[1, 1] =
				(jacobianMatrix[0, 0] * jacobianMatrix[2, 2] - jacobianMatrix[0, 2] * jacobianMatrix[2, 0]) / jacdet;
			inverseJacobian[1, 2] =
				(jacobianMatrix[0, 2] * jacobianMatrix[1, 0] - jacobianMatrix[0, 0] * jacobianMatrix[1, 2]) / jacdet;

			inverseJacobian[2, 0] =
				(jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[1, 1] * jacobianMatrix[2, 0]) / jacdet;
			inverseJacobian[2, 1] =
				(jacobianMatrix[0, 1] * jacobianMatrix[2, 0] - jacobianMatrix[0, 0] * jacobianMatrix[2, 1]) / jacdet;
			inverseJacobian[2, 2] =
				(jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0]) / jacdet;

			return inverseJacobian;
		}

		private static double CalculateJacobianDeterminant(double[,] jacobianMatrix)
		{
			return jacobianMatrix[0, 0] *
			       (jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[2, 1] * jacobianMatrix[1, 2])
			       - jacobianMatrix[0, 1] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 2] -
			                                 jacobianMatrix[2, 0] * jacobianMatrix[1, 2])
			       + jacobianMatrix[0, 2] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 1] -
			                                 jacobianMatrix[2, 0] * jacobianMatrix[1, 1]);
		}

		private static double[,] CalculateJacobian(ControlPoint[] elementControlPoints, NURBS3D nurbs, int j)
		{
			var jacobianMatrix = new double[3, 3];

			for (int k = 0; k < elementControlPoints.Length; k++)
			{
				jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
				jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[2, 0] += nurbs.NurbsDerivativeValuesZeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[2, 1] += nurbs.NurbsDerivativeValuesZeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[2, 2] += nurbs.NurbsDerivativeValuesZeta[k, j] * elementControlPoints[k].Z;
			}

			return jacobianMatrix;
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta,
				element.Patch.DegreeZeta, element.Knots);
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
			NeumannBoundaryCondition neumann)
		{
			throw new NotSupportedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
			NeumannBoundaryCondition neumann)
		{
			throw new NotSupportedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
			PressureBoundaryCondition pressure)
		{
			throw new NotSupportedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
			PressureBoundaryCondition pressure)
		{
			throw new NotSupportedException();
		}

		public double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements)
		{
			var nurbsElement = (NURBSElement3D) element;
			var knotParametricCoordinatesKsi = new Vector(new double[] {element.Knots[0].Ksi, element.Knots[4].Ksi});
			var knotParametricCoordinatesHeta = new Vector(new double[] {element.Knots[0].Heta, element.Knots[2].Heta});
			var knotParametricCoordinatesΖeta = new Vector(new double[] {element.Knots[0].Zeta, element.Knots[1].Zeta});
			NURBS3D nurbs = new NURBS3D(nurbsElement, nurbsElement.ControlPoints, knotParametricCoordinatesKsi,
				knotParametricCoordinatesHeta, knotParametricCoordinatesΖeta);
			var knotDisplacements = new double[8, 3];
			var paraviewKnotRenumbering = new int[] {0, 4, 2, 6, 1, 5, 3, 7};
			for (int j = 0; j < element.Knots.Count; j++)
			{
				for (int i = 0; i < element.ControlPoints.Count; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] +=
						nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
					knotDisplacements[paraviewKnotRenumbering[j], 1] +=
						nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
					knotDisplacements[paraviewKnotRenumbering[j], 2] +=
						nurbs.NurbsValues[i, j] * localDisplacements[i, 2];
				}
			}

			return knotDisplacements;
		}
	}
}