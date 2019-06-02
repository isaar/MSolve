using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Elements
{
    public class TSplineKirchhoffLoveShellSectionElement: Element, IStructuralIsogeometricElement
	{
	    public Matrix ExtractionOperator { get; set; } 
		public int DegreeKsi { get; set; }
		public int DegreeHeta { get; set; }
		protected readonly static IDofType[] controlPointDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		protected IDofType[][] dofTypes;
		protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private DynamicMaterial dynamicProperties;
		private IReadOnlyList<IShellSectionMaterial> materialsAtGaussPoints;
        public CellType CellType { get; } = CellType.Unknown;

        public IElementDofEnumerator DofEnumerator
		{
			get
			{
				return dofEnumerator;
			}

			set
			{
				this.dofEnumerator = value;
			}
		}

		public ElementDimensions ElementDimensions
		{
			get
			{
				return ElementDimensions.ThreeD;
			}
		}

		public bool MaterialModified => throw new NotImplementedException();

		public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
		{
			var shellElement = (TSplineKirchhoffLoveShellElement)element;
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(shellElement);
			var ElementNodalForces = new double[shellElement.ControlPointsDictionary.Count * 3];
			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);
			
			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

				var hessianMatrix = CalculateHessian(shellElement, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = surfaceBasisVector1 .CrossProduct( surfaceBasisVector2);
				var J1 = surfaceBasisVector3.Norm2();
				surfaceBasisVector3.ScaleIntoThis(1 / J1);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2, surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, shellElement);

				
				var MembraneForces = materialsAtGaussPoints[j].MembraneForces;
				var BendingMoments = materialsAtGaussPoints[j].Moments;


                var Fmembrane = Bmembrane.Multiply(MembraneForces, true);
				Fmembrane.Scale(J1 * gaussPoints[j].WeightFactor);

				var Fbending = Bbending.Multiply(BendingMoments, true);
				Fbending.Scale(J1 * gaussPoints[j].WeightFactor);

				ElementNodalForces.AddIntoThis(Fmembrane);
				ElementNodalForces.AddIntoThis(Fbending);

			}

			return ElementNodalForces;
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
			var shellElement = (TSplineKirchhoffLoveShellElement)element;
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(shellElement);
			//Matrix stiffnessMatrixElement = Matrix.CreateZero(shellElement.ControlPointsDictionary.Count * 3, shellElement.ControlPointsDictionary.Count * 3);

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

				var hessianMatrix = CalculateHessian(shellElement, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = surfaceBasisVector1 .CrossProduct( surfaceBasisVector2);
				var J1 = surfaceBasisVector3.Norm2();
				surfaceBasisVector3.ScaleIntoThis(1 / J1);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12, shellElement);

                var membraneStrain = Bmembrane.Multiply(localDisplacements, true);
				var bendingStrain = Bbending.Multiply(localDisplacements, true);

				materialsAtGaussPoints[j].UpdateMaterial(membraneStrain, bendingStrain);


			}
			return new Tuple<double[], double[]>(new double[0], new double[0]);

		}

		public void ClearMaterialState()
		{
			throw new NotImplementedException();
		}

		public IMatrix DampingMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			var nurbsElement = (TSplineKirchhoffLoveShellElement)element;
			dofTypes = new IDofType[nurbsElement.ControlPoints.Count][];
			for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
			{
				dofTypes[i] = controlPointDOFTypes;
			}
			return dofTypes;
		}

		public IMatrix MassMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public IMatrix StiffnessMatrix(IElement element)
		{
			var shellElement = (TSplineKirchhoffLoveShellElement)element;
            IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(shellElement);
			Matrix stiffnessMatrixElement = Matrix.CreateZero(shellElement.ControlPointsDictionary.Count * 3, shellElement.ControlPointsDictionary.Count * 3);

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

				var hessianMatrix = CalculateHessian(shellElement, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
				var J1 = surfaceBasisVector3.Norm2();
				surfaceBasisVector3.ScaleIntoThis(1 / J1);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);
				
				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);

				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2, surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, shellElement);

				

				var Kmembrane = Bmembrane.Transpose()*(Matrix)materialsAtGaussPoints[j].MembraneConstitutiveMatrix  * Bmembrane  * J1 *
								gaussPoints[j].WeightFactor;

				
				var Kbending = Bbending.Transpose() * (Matrix)materialsAtGaussPoints[j].BendingConstitutiveMatrix * Bbending  * J1 *
							   gaussPoints[j].WeightFactor;

				var KMembraneBending = Bmembrane.Transpose() * (Matrix)materialsAtGaussPoints[j].CouplingConstitutiveMatrix * Bbending * J1 *
				               gaussPoints[j].WeightFactor;

				var KBendingMembrane = Bbending.Transpose() * (Matrix)materialsAtGaussPoints[j].CouplingConstitutiveMatrix * Bmembrane * J1 *
				                       gaussPoints[j].WeightFactor;


				stiffnessMatrixElement.AddIntoThis(Kmembrane);
				stiffnessMatrixElement.AddIntoThis(Kbending);
				stiffnessMatrixElement.AddIntoThis(KMembraneBending);
				stiffnessMatrixElement.AddIntoThis(KBendingMembrane);
			}
			return stiffnessMatrixElement;
		}

		private Matrix CalculateConstitutiveMatrix(TSplineKirchhoffLoveShellElement element, Vector surfaceBasisVector1, Vector surfaceBasisVector2)
		{
            var auxMatrix1 = Matrix.CreateZero(2, 2);
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
			auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
			auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            (Matrix inverse, double det) = auxMatrix1.InvertAndDeterminant();

			var material = ((IContinuumMaterial2D)element.Patch.Material);
			var constitutiveMatrix = Matrix.CreateFromArray(new double[3, 3]
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

		private Matrix CalculateBendingDeformationMatrix(Vector surfaceBasisVector3, ShapeTSplines2DFromBezierExtraction tsplines, int j,
			Vector surfaceBasisVector2, Vector surfaceBasisVectorDerivative1, Vector surfaceBasisVector1, double J1,
			Vector surfaceBasisVectorDerivative2, Vector surfaceBasisVectorDerivative12, TSplineKirchhoffLoveShellElement element)
		{
			Matrix Bbending = Matrix.CreateZero(3, element.ControlPoints.Count * 3);
			for (int column = 0; column < element.ControlPoints.Count * 3; column+=3)
			{
				#region BI1

				var BI1 = surfaceBasisVector3.CrossProduct(surfaceBasisVector3);
				BI1.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				var auxVector = surfaceBasisVector2.CrossProduct(surfaceBasisVector3);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI1.AddIntoThis(auxVector);
				BI1.ScaleIntoThis(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative1));
				auxVector = surfaceBasisVector1.CrossProduct(surfaceBasisVectorDerivative1);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI1.AddIntoThis(auxVector);
				BI1.ScaleIntoThis(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.ScaleIntoThis(-tsplines.TSplineSecondDerivativesValueKsi[column / 3, j]);
				BI1.AddIntoThis(auxVector);

				#endregion

				#region BI2

				Vector BI2 = surfaceBasisVector3.CrossProduct(surfaceBasisVector3);
				BI2.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				auxVector = surfaceBasisVector2.CrossProduct(surfaceBasisVector3);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				BI2.ScaleIntoThis(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative2));
				auxVector = surfaceBasisVector1.CrossProduct(surfaceBasisVectorDerivative2);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				auxVector = surfaceBasisVectorDerivative2.CrossProduct(surfaceBasisVector2);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				BI2.ScaleIntoThis(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.ScaleIntoThis(-tsplines.TSplineSecondDerivativesValueHeta[column / 3, j]);
				BI2.AddIntoThis(auxVector);

				#endregion

				#region BI3

				Vector BI3 = surfaceBasisVector3.CrossProduct(surfaceBasisVector3);
				BI3.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				auxVector = surfaceBasisVector2.CrossProduct(surfaceBasisVector3);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				BI3.ScaleIntoThis(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative12));
				auxVector = surfaceBasisVector1.CrossProduct(surfaceBasisVectorDerivative12);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				auxVector = surfaceBasisVectorDerivative2.CrossProduct(surfaceBasisVector2);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				BI3.ScaleIntoThis(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.ScaleIntoThis(-tsplines.TSplineSecondDerivativesValueKsiHeta[column / 3, j]);
				BI3.AddIntoThis(auxVector);

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

		private Matrix CalculateMembraneDeformationMatrix(ShapeTSplines2DFromBezierExtraction tsplines, int j, Vector surfaceBasisVector1,
			Vector surfaceBasisVector2, TSplineKirchhoffLoveShellElement element)
		{
			Matrix dRIa = Matrix.CreateZero(3, element.ControlPoints.Count * 3);
			for (int i = 0; i < element.ControlPoints.Count; i++)
			{
				for (int m = 0; m < 3; m++)
				{
					dRIa[m, i] = tsplines.TSplineDerivativeValuesHeta[i, j] * surfaceBasisVector1[m] +
								 tsplines.TSplineDerivativeValuesKsi[i, j] * surfaceBasisVector2[m];
				}
			}

			Matrix Bmembrane = Matrix.CreateZero(3, element.ControlPoints.Count * 3);
			for (int column = 0; column < element.ControlPoints.Count * 3; column+=3)
			{
				Bmembrane[0, column] = tsplines.TSplineDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[0];
				Bmembrane[0, column + 1] = tsplines.TSplineDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[1];
				Bmembrane[0, column + 2] = tsplines.TSplineDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[2];

				Bmembrane[1, column] = tsplines.TSplineDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[0];
				Bmembrane[1, column + 1] = tsplines.TSplineDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[1];
				Bmembrane[1, column + 2] = tsplines.TSplineDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[2];

				Bmembrane[2, column] = dRIa[0, column / 3];
				Bmembrane[2, column + 1] = dRIa[1, column / 3];
				Bmembrane[2, column + 2] = dRIa[2, column / 3];
			}

			return Bmembrane;
		}

		private static Vector CalculateSurfaceBasisVector1(Matrix Matrix, int row)
		{
			Vector surfaceBasisVector1 = Vector.CreateZero(3);
			surfaceBasisVector1[0] = Matrix[row, 0];
			surfaceBasisVector1[1] = Matrix[row, 1];
			surfaceBasisVector1[2] = Matrix[row, 2];
			return surfaceBasisVector1;
		}

		private static Matrix CalculateHessian(TSplineKirchhoffLoveShellElement shellElement, ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			Matrix hessianMatrix = Matrix.CreateZero(3, 3);
			for (int k = 0; k < shellElement.ControlPoints.Count; k++)
			{
				hessianMatrix[0, 0] += tsplines.TSplineSecondDerivativesValueKsi[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[0, 1] += tsplines.TSplineSecondDerivativesValueKsi[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[0, 2] += tsplines.TSplineSecondDerivativesValueKsi[k, j] * shellElement.ControlPoints[k].Z;
				hessianMatrix[1, 0] += tsplines.TSplineSecondDerivativesValueHeta[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[1, 1] += tsplines.TSplineSecondDerivativesValueHeta[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[1, 2] += tsplines.TSplineSecondDerivativesValueHeta[k, j] * shellElement.ControlPoints[k].Z;
				hessianMatrix[2, 0] += tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[2, 1] += tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[2, 2] += tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * shellElement.ControlPoints[k].Z;
			}

			return hessianMatrix;
		}

		private static Matrix CalculateJacobian(TSplineKirchhoffLoveShellElement shellElement, ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
			for (int k = 0; k < shellElement.ControlPoints.Count; k++)
			{
				jacobianMatrix[0, 0] += tsplines.TSplineDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].X;
				jacobianMatrix[0, 1] += tsplines.TSplineDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].Y;
				jacobianMatrix[0, 2] += tsplines.TSplineDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].Z;
				jacobianMatrix[1, 0] += tsplines.TSplineDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].X;
				jacobianMatrix[1, 1] += tsplines.TSplineDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].Y;
				jacobianMatrix[1, 2] += tsplines.TSplineDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].Z;
			}

			return jacobianMatrix;
		}

        private IList<GaussLegendrePoint3D> CreateElementGaussPoints(TSplineKirchhoffLoveShellElement element)
        {
            GaussQuadrature gauss = new GaussQuadrature();
            return gauss.CalculateElementGaussPoints(element.DegreeKsi, element.DegreeHeta, new List<Knot>
                {
                    new Knot(){ID=0,Ksi=-1,Heta = -1,Zeta = 0},
                    new Knot(){ID=1,Ksi=-1,Heta = 1,Zeta = 0},
                    new Knot(){ID=2,Ksi=1,Heta = -1,Zeta = 0},
                    new Knot(){ID=3,Ksi=1,Heta = 1,Zeta = 0}
                });
		}

		public double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements)
		{
			throw new NotImplementedException();
		}
	}
}
