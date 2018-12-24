using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.IGA.Elements
{
    public class TSplineKirchhoffLoveShellElementMaterial: Element, IStructuralIsogeometricElement
	{
	    public Matrix2D ExtractionOperator { get; set; } 
		public int DegreeKsi { get; set; }
		public int DegreeHeta { get; set; }
		protected readonly static DOFType[] controlPointDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
		protected DOFType[][] dofTypes;
		protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
		private DynamicMaterial dynamicProperties;
		private IReadOnlyList<IShellMaterial> materialsAtGaussPoints;

		private Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>
			materialsAtThicknessGP =
				new Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>();
		

		public TSplineKirchhoffLoveShellElementMaterial(int id, Patch patch, int degreeKsi, int degreeHeta,
			double thickness, Matrix2D extractionOperator, IShellMaterial shellMaterial)
		{
			this.ID = id;
			this.Patch = patch;
			this.DegreeKsi = degreeKsi;
			this.DegreeHeta = degreeHeta;
			this.Thickness = thickness;
			this.ExtractionOperator = extractionOperator;

			CreateElementGaussPoints(this);
			foreach (var medianSurfaceGP in thicknessIntegrationPoints.Keys)
			{
				materialsAtThicknessGP.Add(medianSurfaceGP, new Dictionary<GaussLegendrePoint3D, IShellMaterial>());
				foreach (var point in thicknessIntegrationPoints[medianSurfaceGP])
				{
					materialsAtThicknessGP[medianSurfaceGP].Add(point, shellMaterial.Clone());
				}
			}
		}

		public IElementDOFEnumerator DOFEnumerator
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

		public double Thickness { get; set; }

		public bool MaterialModified => throw new NotImplementedException();

		public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
		{
			var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(shellElement);
			Vector ElementNodalForces = new Vector(shellElement.ControlPointsDictionary.Count * 3);
			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);


			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

				var hessianMatrix = CalculateHessian(shellElement, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = surfaceBasisVector1 ^ surfaceBasisVector2;
				var J1 = surfaceBasisVector3.Norm;
				surfaceBasisVector3.Multiply(1 / J1);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2, surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, shellElement);				

				var (MembraneForces, BendingMoments) =
					IntegratedStressesOverThickness(gaussPoints, j);

				var Fmembrane = Bmembrane.Transpose() * MembraneForces;
				Fmembrane.Scale(J1 * gaussPoints[j].WeightFactor);

				var Fbending = Bbending.Transpose() * BendingMoments;
				Fbending.Scale(J1 *gaussPoints[j].WeightFactor);
				               				
				ElementNodalForces.Add(Fmembrane);
				ElementNodalForces.Add(Fbending);
				
			}

			return ElementNodalForces.Data;
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
			var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			//Matrix2D stiffnessMatrixElement = new Matrix2D(shellElement.ControlPointsDictionary.Count * 3, shellElement.ControlPointsDictionary.Count * 3);

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);

			for (int j = 0; j < materialsAtThicknessGP.Keys.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

				var hessianMatrix = CalculateHessian(shellElement, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = surfaceBasisVector1 ^ surfaceBasisVector2;
				var J1 = surfaceBasisVector3.Norm;
				surfaceBasisVector3.Multiply(1 / J1);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12, shellElement);

				var membraneStrain = Bmembrane.Transpose() * new Vector(localDisplacements);
				var bendingStrain = Bbending.Transpose() * (new Vector(localDisplacements));
				
				foreach (var keyValuePair in materialsAtThicknessGP[materialsAtThicknessGP.Keys.ToList()[j]])
				{
					var thicknessPoint = keyValuePair.Key;
					var material = keyValuePair.Value;

					var bendingStrainContribution= new Vector(bendingStrain.Length);
					bendingStrain.CopyTo(bendingStrainContribution.Data,0);
					bendingStrainContribution.Scale(-thicknessPoint.Zeta);
					var gpStrain = membraneStrain + bendingStrainContribution;
					material.UpdateMaterial(gpStrain);
				}


			}
			return new Tuple<double[], double[]>(new double[0], new double[0]);

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
			var nurbsElement = (TSplineKirchhoffLoveShellElementMaterial)element;
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
			var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			Matrix2D stiffnessMatrixElement = new Matrix2D(shellElement.ControlPointsDictionary.Count * 3, shellElement.ControlPointsDictionary.Count * 3);
			var gaussPoints = materialsAtThicknessGP.Keys.ToList();
			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

				var hessianMatrix = CalculateHessian(shellElement, tsplines, j);
				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = surfaceBasisVector1^surfaceBasisVector2;

				foreach (var integrationPointMaterial in materialsAtThicknessGP[gaussPoints[j]].Values)
				{
					integrationPointMaterial.TangentVectorV1 = surfaceBasisVector1.Data;
					integrationPointMaterial.TangentVectorV2 = surfaceBasisVector2.Data;
					integrationPointMaterial.NormalVectorV3 = surfaceBasisVector3.Data;
				}
				var J1 = surfaceBasisVector3.Norm;
				surfaceBasisVector3.Multiply(1 / J1);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12, shellElement);

				var (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix) =
					IntegratedConstitutiveOverThickness(gaussPoints, j);

				var Kmembrane = Bmembrane.Transpose() * MembraneConstitutiveMatrix * Bmembrane  * J1 *
				                gaussPoints[j].WeightFactor;
				var Kbending = Bbending.Transpose() * BendingConstitutiveMatrix * Bbending  * J1 *
							   gaussPoints[j].WeightFactor;

				var KMembraneBending = Bmembrane.Transpose() * CouplingConstitutiveMatrix * Bbending * J1 *
				                       gaussPoints[j].WeightFactor;

				var KBendingMembrane = Bbending.Transpose() * CouplingConstitutiveMatrix * Bmembrane * J1 *
				                       gaussPoints[j].WeightFactor;


				stiffnessMatrixElement.Add(Kmembrane);
				stiffnessMatrixElement.Add(Kbending);
				stiffnessMatrixElement.Add(KMembraneBending);
				stiffnessMatrixElement.Add(KBendingMembrane);
			}
			return stiffnessMatrixElement;
		}

		public (Matrix2D MembraneConstitutiveMatrix, Matrix2D BendingConstitutiveMatrix, Matrix2D CouplingConstitutiveMatrix) IntegratedConstitutiveOverThickness(IList<GaussLegendrePoint3D> gaussPoints, int j)
		{
			var MembraneConstitutiveMatrix = new Matrix2D(3, 3);
			var BendingConstitutiveMatrix = new Matrix2D(3, 3);
			var CouplingConstitutiveMatrix = new Matrix2D(3, 3);

			foreach (var keyValuePair in materialsAtThicknessGP[gaussPoints[j]])
			{
				var thicknessPoint = keyValuePair.Key;
				var material = keyValuePair.Value;
				var constitutiveMatrixM= CopyConstitutiveMatrix((Matrix2D)material.ConstitutiveMatrix);
				MembraneConstitutiveMatrix.Add(constitutiveMatrixM * thicknessPoint.WeightFactor);
				var constitutiveMatrixB = CopyConstitutiveMatrix((Matrix2D)material.ConstitutiveMatrix);
				BendingConstitutiveMatrix.Add(constitutiveMatrixB * thicknessPoint.WeightFactor *
											  Math.Pow(thicknessPoint.Zeta, 2));
				var constitutiveMatrixC = CopyConstitutiveMatrix((Matrix2D)material.ConstitutiveMatrix);
				CouplingConstitutiveMatrix.Add(constitutiveMatrixC * thicknessPoint.WeightFactor *
											   thicknessPoint.Zeta);
			}

			return (MembraneConstitutiveMatrix, BendingConstitutiveMatrix,CouplingConstitutiveMatrix);
		}

		private Matrix2D CopyConstitutiveMatrix(Matrix2D f)
		{
			double [,] g= new double[f.Rows,f.Columns];
			Array.Copy(f.Data, 0, g, 0, f.Data.Length);
			return new Matrix2D(g);
		}

		public (Vector MembraneForces, Vector BendingMoments ) IntegratedStressesOverThickness(
			IList<GaussLegendrePoint3D> gaussPoints, int j)
		{
			var MembraneForces = new Vector(3);
			var BendingMoments = new Vector(3);

			foreach (var keyValuePair in materialsAtThicknessGP[gaussPoints[j]])
			{
				var thicknessPoint = keyValuePair.Key;
				var material = keyValuePair.Value;
				var gpForces = new Vector(material.Stresses);
				gpForces.Scale(thicknessPoint.WeightFactor * (Thickness / 2));
				MembraneForces.Add(gpForces);

				var gpMoments = new Vector(material.Stresses);
				gpMoments.Scale(thicknessPoint.WeightFactor *
				                Math.Pow(thicknessPoint.Zeta, 2) * (Thickness / 2));
				BendingMoments.Add(gpMoments);

			}

			return (MembraneForces, BendingMoments);
		}

		
		private Matrix2D CalculateBendingDeformationMatrix(Vector surfaceBasisVector3, ShapeTSplines2DFromBezierExtraction tsplines, int j,
			Vector surfaceBasisVector2, Vector surfaceBasisVectorDerivative1, Vector surfaceBasisVector1, double J1,
			Vector surfaceBasisVectorDerivative2, Vector surfaceBasisVectorDerivative12, TSplineKirchhoffLoveShellElementMaterial element)
		{
			Matrix2D Bbending = new Matrix2D(3, element.ControlPoints.Count * 3);
			for (int column = 0; column < element.ControlPoints.Count * 3; column+=3)
			{
				#region BI1

				var BI1 = surfaceBasisVector3^surfaceBasisVector3;
				BI1.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				var auxVector = surfaceBasisVector2^surfaceBasisVector3;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI1.Add(auxVector);
				BI1.Multiply(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative1));
				auxVector = surfaceBasisVector1^surfaceBasisVectorDerivative1;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI1.Add(auxVector);
				BI1.Multiply(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.Multiply(-tsplines.TSplineSecondDerivativesValueKsi[column / 3, j]);
				BI1.Add(auxVector);

				#endregion

				#region BI2

				Vector BI2 = surfaceBasisVector3^surfaceBasisVector3;
				BI2.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				auxVector = surfaceBasisVector2^surfaceBasisVector3;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI2.Add(auxVector);
				BI2.Multiply(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative2));
				auxVector = surfaceBasisVector1^surfaceBasisVectorDerivative2;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI2.Add(auxVector);
				auxVector = surfaceBasisVectorDerivative2^surfaceBasisVector2;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI2.Add(auxVector);
				BI2.Multiply(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.Multiply(-tsplines.TSplineSecondDerivativesValueHeta[column / 3, j]);
				BI2.Add(auxVector);

				#endregion

				#region BI3

				Vector BI3 = surfaceBasisVector3^surfaceBasisVector3;
				BI3.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				auxVector = surfaceBasisVector2^surfaceBasisVector3;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI3.Add(auxVector);
				BI3.Multiply(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative12));
				auxVector = surfaceBasisVector1^surfaceBasisVectorDerivative12;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI3.Add(auxVector);
				auxVector = surfaceBasisVectorDerivative2^surfaceBasisVector2;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI3.Add(auxVector);
				BI3.Multiply(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.Multiply(-tsplines.TSplineSecondDerivativesValueKsiHeta[column / 3, j]);
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

		private Matrix2D CalculateMembraneDeformationMatrix(ShapeTSplines2DFromBezierExtraction tsplines, int j, Vector surfaceBasisVector1,
			Vector surfaceBasisVector2, TSplineKirchhoffLoveShellElementMaterial element)
		{
			Matrix2D dRIa = new Matrix2D(3, element.ControlPoints.Count * 3);
			for (int i = 0; i < element.ControlPoints.Count; i++)
			{
				for (int m = 0; m < 3; m++)
				{
					dRIa[m, i] = tsplines.TSplineDerivativeValuesHeta[i, j] * surfaceBasisVector1[m] +
								 tsplines.TSplineDerivativeValuesKsi[i, j] * surfaceBasisVector2[m];
				}
			}

			Matrix2D Bmembrane = new Matrix2D(3, element.ControlPoints.Count * 3);
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

		private static Vector CalculateSurfaceBasisVector1(Matrix2D Matrix, int row)
		{
			Vector surfaceBasisVector1 = new Vector(3);
			surfaceBasisVector1[0] = Matrix[row, 0];
			surfaceBasisVector1[1] = Matrix[row, 1];
			surfaceBasisVector1[2] = Matrix[row, 2];
			return surfaceBasisVector1;
		}

		private static Matrix2D CalculateHessian(TSplineKirchhoffLoveShellElementMaterial shellElement, ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			Matrix2D hessianMatrix = new Matrix2D(3, 3);
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

		private static Matrix2D CalculateJacobian(TSplineKirchhoffLoveShellElementMaterial shellElement, ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			Matrix2D jacobianMatrix = new Matrix2D(2, 3);
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


		const int thicknessIntegrationDegree = 2;
		Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>> thicknessIntegrationPoints= new Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>>();
		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(TSplineKirchhoffLoveShellElementMaterial element)
        {
            GaussQuadrature gauss = new GaussQuadrature();
	        var medianSurfaceGP= gauss.CalculateElementGaussPoints(element.DegreeKsi, element.DegreeHeta, new List<Knot>
	        {
		        new Knot() {ID = 0, Ksi = -1, Heta = -1, Zeta = 0},
		        new Knot() {ID = 1, Ksi = -1, Heta = 1, Zeta = 0},
		        new Knot() {ID = 2, Ksi = 1, Heta = -1, Zeta = 0},
		        new Knot() {ID = 3, Ksi = 1, Heta = 1, Zeta = 0},
	        });

	        foreach (var point in medianSurfaceGP)
	        {
				thicknessIntegrationPoints.Add(point, gauss.CalculateElementGaussPoints(thicknessIntegrationDegree,
					new List<Knot>
					{
						new Knot() {ID = 0, Ksi = -Thickness / 2, Heta = point.Heta },
						new Knot() {ID = 1, Ksi = Thickness / 2, Heta = point.Heta},
					}).Select(gp => new GaussLegendrePoint3D(point.Ksi, point.Heta, gp.Ksi, null, gp.WeightFactor)).ToList());
	        }

	        //new Knot() { ID = 0, Ksi = -1, Heta = -1, Zeta = -Thickness / 2 },
	        //new Knot() { ID = 1, Ksi = -1, Heta = -1, Zeta = Thickness / 2 },
	        //new Knot() { ID = 2, Ksi = -1, Heta = 1, Zeta = -Thickness / 2 },
	        //new Knot() { ID = 3, Ksi = -1, Heta = 1, Zeta = Thickness / 2 },
	        //new Knot() { ID = 4, Ksi = 1, Heta = -1, Zeta = -Thickness / 2 },
	        //new Knot() { ID = 5, Ksi = 1, Heta = -1, Zeta = Thickness / 2 },
	        //new Knot() { ID = 6, Ksi = 1, Heta = 1, Zeta = -Thickness / 2 },
	        //new Knot() { ID = 7, Ksi = 1, Heta = 1, Zeta = Thickness / 2 }


	        return medianSurfaceGP;
        }
	}
}
