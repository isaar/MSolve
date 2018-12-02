using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.IGA.Elements
{
	public class NURBSKirchhoffLoveShellElement : Element, IStructuralIsogeometricElement
	{
		protected readonly static DOFType[] controlPointDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
		protected DOFType[][] dofTypes;
		protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();





		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IElementDOFEnumerator DOFEnumerator
		{
			get => dofEnumerator;
			set => this.dofEnumerator = value;
		}

		public bool MaterialModified => throw new NotImplementedException();

		public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure)
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		public void ClearMaterialState()
		{
			throw new NotImplementedException();
		}

		public IMatrix2D DampingMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
		{
			var nurbsElement = (NURBSKirchhoffLoveShellElement)element;
			dofTypes = new DOFType[nurbsElement.ControlPoints.Count][];
			for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
			{
				dofTypes[i] = controlPointDOFTypes;
			}
			return dofTypes;
		}

		public IMatrix2D MassMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public IMatrix2D StiffnessMatrix(IElement element)
		{
			var shellElement = (NURBSKirchhoffLoveShellElement) element;
			
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(shellElement);
			Matrix2D stiffnessMatrixElement = new Matrix2D(shellElement.ControlPointsDictionary.Count * 3, shellElement.ControlPointsDictionary.Count * 3);

			NURBS2D nurbs = new NURBS2D(shellElement, shellElement.ControlPoints);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, nurbs, j);

				var hessianMatrix = CalculateHessian(shellElement, nurbs, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix,0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = surfaceBasisVector1^surfaceBasisVector2;
				var J1 = surfaceBasisVector3.Norm;
				surfaceBasisVector3.Multiply(1 / J1);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				Matrix2D constitutiveMatrix = CalculateConstitutiveMatrix(shellElement, surfaceBasisVector1, surfaceBasisVector2);

				 var Bmembrane = CalculateMembraneDeformationMatrix(nurbs, j, surfaceBasisVector1, surfaceBasisVector2,shellElement);

				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, nurbs, j, surfaceBasisVector2, surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, shellElement);

				double membraneStiffness = ((IIsotropicContinuumMaterial2D)shellElement.Patch.Material).YoungModulus * shellElement.Patch.Thickness /
				                           (1 - Math.Pow(((IIsotropicContinuumMaterial2D)shellElement.Patch.Material).PoissonRatio, 2));

				var Kmembrane = Bmembrane.Transpose() * constitutiveMatrix * Bmembrane * membraneStiffness * J1 *
				                gaussPoints[j].WeightFactor;

				double bendingStiffness = ((IIsotropicContinuumMaterial2D)shellElement.Patch.Material).YoungModulus * Math.Pow(shellElement.Patch.Thickness, 3) /
				                          12 / (1 - Math.Pow(((IIsotropicContinuumMaterial2D)shellElement.Patch.Material).PoissonRatio, 2));

				var Kbending = Bbending.Transpose() * constitutiveMatrix * Bbending * bendingStiffness * J1 *
				               gaussPoints[j].WeightFactor;


				stiffnessMatrixElement.Add(Kmembrane);
				stiffnessMatrixElement.Add(Kbending);
			}
			return stiffnessMatrixElement;
		}

		private Matrix2D CalculateConstitutiveMatrix(NURBSKirchhoffLoveShellElement element,Vector surfaceBasisVector1, Vector surfaceBasisVector2)
		{
            var auxMatrix1 = new Matrix2D(2, 2);
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
			auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            (Matrix2D inverse, double det) = auxMatrix1.Invert2x2AndDeterminant();

			var material =((IContinuumMaterial2D)element.Patch.Material);
			var constitutiveMatrix = new Matrix2D(new double[3, 3]
			{
				{
                    inverse[0,0]*inverse[0,0],
					material.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-material.PoissonRatio)*inverse[1,0]*inverse[1,0],
					inverse[0,0]*inverse[1,0]
				},
				{
					material.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-material.PoissonRatio)*inverse[1,0]*inverse[1,0],
					inverse[1,1]*inverse[1,1],
					inverse[1,1]*inverse[1,0]
				},
				{
					inverse[0,0]*inverse[1,0],
					inverse[1,1]*inverse[1,0],
					0.5*(1-material.PoissonRatio)*inverse[0,0]*inverse[1,1]+(1+material.PoissonRatio)*inverse[1,0]*inverse[1,0]
				},
			});
			return constitutiveMatrix;
		}

		private Matrix2D CalculateBendingDeformationMatrix(Vector surfaceBasisVector3, NURBS2D nurbs, int j,
			Vector surfaceBasisVector2, Vector surfaceBasisVectorDerivative1, Vector surfaceBasisVector1, double J1,
			Vector surfaceBasisVectorDerivative2, Vector surfaceBasisVectorDerivative12, NURBSKirchhoffLoveShellElement element)
		{
			Matrix2D Bbending = new Matrix2D(3, element.ControlPoints.Count * 3);
			for (int column = 0; column < element.ControlPoints.Count * 3; column+=3)
			{
				#region BI1

				var BI1 = surfaceBasisVector3^surfaceBasisVector3;
				BI1.Multiply(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
				var auxVector = surfaceBasisVector2^surfaceBasisVector3;
				auxVector.Multiply(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
				BI1.Add(auxVector);
				BI1.Multiply(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative1));
				auxVector = surfaceBasisVector1^surfaceBasisVectorDerivative1;
				auxVector.Multiply(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
				BI1.Add(auxVector);
				BI1.Multiply(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.Multiply(-nurbs.NurbsSecondDerivativeValueKsi[column / 3, j]);
				BI1.Add(auxVector);

				#endregion

				#region BI2

				Vector BI2 = surfaceBasisVector3^surfaceBasisVector3;
				BI2.Multiply(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
				auxVector = surfaceBasisVector2^surfaceBasisVector3;
				auxVector.Multiply(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
				BI2.Add(auxVector);
				BI2.Multiply(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative2));
				auxVector = surfaceBasisVector1^surfaceBasisVectorDerivative2;
				auxVector.Multiply(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
				BI2.Add(auxVector);
				auxVector = surfaceBasisVectorDerivative2^surfaceBasisVector2;
				auxVector.Multiply(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
				BI2.Add(auxVector);
				BI2.Multiply(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.Multiply(-nurbs.NurbsSecondDerivativeValueHeta[column / 3, j]);
				BI2.Add(auxVector);

				#endregion

				#region BI3

				Vector BI3 = surfaceBasisVector3^surfaceBasisVector3;
				BI3.Multiply(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
				auxVector = surfaceBasisVector2^surfaceBasisVector3;
				auxVector.Multiply(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
				BI3.Add(auxVector);
				BI3.Multiply(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative12));
				auxVector = surfaceBasisVector1^surfaceBasisVectorDerivative12;
				auxVector.Multiply(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
				BI3.Add(auxVector);
				auxVector = surfaceBasisVectorDerivative2^surfaceBasisVector2;
				auxVector.Multiply(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
				BI3.Add(auxVector);
				BI3.Multiply(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.Multiply(-nurbs.NurbsSecondDerivativeValueKsiHeta[column / 3, j]);
				BI3.Add(auxVector);

				#endregion

				Bbending[0, column] = BI1[0];
				Bbending[0, column + 1] = BI1[1];
				Bbending[0, column + 2] = BI1[2];

				Bbending[1, column] = BI2[0];
				Bbending[1, column + 1] = BI2[1];
				Bbending[1, column + 2] = BI2[2];

				Bbending[2, column] = 2 * BI3[0];
				Bbending[2, column + 1] = 2 * BI3[1];
				Bbending[2, column + 2] = 2 * BI3[2];
			}

			return Bbending;
		}

		private Matrix2D CalculateMembraneDeformationMatrix(NURBS2D nurbs, int j, Vector surfaceBasisVector1,
			Vector surfaceBasisVector2, NURBSKirchhoffLoveShellElement element)
		{
			Matrix2D dRIa = new Matrix2D(3, element.ControlPoints.Count );
			for (int i = 0; i < element.ControlPoints.Count ; i++)
			{
				for (int m = 0; m < 3; m++)
				{
					dRIa[m, i] = nurbs.NurbsDerivativeValuesHeta[i, j] * surfaceBasisVector1[m] +
					             nurbs.NurbsDerivativeValuesKsi[i, j] * surfaceBasisVector2[m];
				}
			}

			Matrix2D Bmembrane = new Matrix2D(3, element.ControlPoints.Count * 3);
			for (int column = 0; column < element.ControlPoints.Count * 3; column+=3)
			{
				Bmembrane[0, column] = nurbs.NurbsDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[0];
				Bmembrane[0, column + 1] = nurbs.NurbsDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[1];
				Bmembrane[0, column + 2] = nurbs.NurbsDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[2];

				Bmembrane[1, column] = nurbs.NurbsDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[0];
				Bmembrane[1, column + 1] = nurbs.NurbsDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[1];
				Bmembrane[1, column + 2] = nurbs.NurbsDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[2];

				Bmembrane[2, column] = dRIa[0, column / 3];
				Bmembrane[2, column + 1] = dRIa[1, column / 3];
				Bmembrane[2, column + 2] = dRIa[2, column / 3];
			}

			return Bmembrane;
		}

		private static Vector CalculateSurfaceBasisVector1(Matrix2D Matrix,  int row)
		{
			Vector surfaceBasisVector1 = new Vector(3);
			surfaceBasisVector1[0] = Matrix[row, 0];
			surfaceBasisVector1[1] = Matrix[row, 1];
			surfaceBasisVector1[2] = Matrix[row, 2];
			return surfaceBasisVector1;
		}

		private static Matrix2D CalculateHessian(NURBSKirchhoffLoveShellElement shellElement, NURBS2D nurbs, int j)
		{
			Matrix2D hessianMatrix = new Matrix2D(3, 3);
			for (int k = 0; k < shellElement.ControlPoints.Count; k++)
			{
				hessianMatrix[0, 0] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[0, 1] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[0, 2] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * shellElement.ControlPoints[k].Z;
				hessianMatrix[1, 0] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[1, 1] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[1, 2] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * shellElement.ControlPoints[k].Z;
				hessianMatrix[2, 0] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[2, 1] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[2, 2] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * shellElement.ControlPoints[k].Z;
			}

			return hessianMatrix;
		}

		private static Matrix2D CalculateJacobian(NURBSKirchhoffLoveShellElement shellElement, NURBS2D nurbs, int j)
		{
			Matrix2D jacobianMatrix = new Matrix2D(2, 3);
			for (int k = 0; k < shellElement.ControlPoints.Count; k++)
			{
				jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].X;
				jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].Y;
				jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].Z;
				jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].X;
				jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].Y;
				jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].Z;
			}

			return jacobianMatrix;
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Knots);
		}
	}
}
