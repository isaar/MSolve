using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.IGA.Elements
{
	public class TSplineElement2D : Element, IStructuralIsogeometricElement
	{
		protected readonly static DOFType[] controlPointDOFTypes = new DOFType[] { DOFType.X, DOFType.Y };
		protected DOFType[][] dofTypes;
		private IReadOnlyList<IContinuumMaterial2D> materialsAtGaussPoints;
		private DynamicMaterial dynamicProperties;

		protected IElementDofEnumerator_v2 dofEnumerator = new GenericDofEnumerator_v2();
		public IElementDofEnumerator_v2 DofEnumerator
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
				return ElementDimensions.TwoD;
			}
		}

		public bool MaterialModified => throw new NotImplementedException();

		public Matrix2D ExtractionOperator { get; set; }

		public int DegreeKsi { get; set; }
		public int DegreeHeta { get; set; }
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

		public IMatrix DampingMatrix(IElement_v2 element)
		{
			throw new NotImplementedException();
		}

		public IList<IList<DOFType>> GetElementDOFTypes(IElement_v2 element)
		{
			var nurbsElement = (TSplineElement2D)element;
			dofTypes = new DOFType[nurbsElement.ControlPoints.Count][];
			for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
			{
				dofTypes[i] = controlPointDOFTypes;
			}
			return dofTypes;
		}

		public IMatrix MassMatrix(IElement_v2 element)
		{
			throw new NotImplementedException();
		}

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public IMatrix StiffnessMatrix(IElement_v2 element)
		{
			var tsplineElement = (TSplineElement2D)element;
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(tsplineElement);
			Matrix2D stiffnessMatrixElement = new Matrix2D(tsplineElement.ControlPointsDictionary.Count * 2, tsplineElement.ControlPointsDictionary.Count * 2);

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(tsplineElement, tsplineElement.ControlPoints);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobianMatrix(tsplineElement, tsplines, j);

				var jacdet = CalculateJacobianDeterminant(jacobianMatrix);

				var B1 = CalculateDeformationMatrix1(jacobianMatrix, jacdet);

				var B2 = CalculateDeformationMatrix2(tsplineElement, tsplines, j);

				Matrix2D B = B1 * B2;
				Matrix2D ElasticityMatrix = materialsAtGaussPoints[j].ConstitutiveMatrix;
				Matrix2D stiffnessMatrixGaussPoint = B.Transpose() * ElasticityMatrix;
				stiffnessMatrixGaussPoint = stiffnessMatrixGaussPoint * B;
				stiffnessMatrixGaussPoint = stiffnessMatrixGaussPoint * (jacdet * gaussPoints[j].WeightFactor * tsplineElement.Patch.Thickness);

				for (int m = 0; m < tsplineElement.ControlPoints.Count * 2; m++)
				{
					for (int n = 0; n < tsplineElement.ControlPoints.Count * 2; n++)
					{
						stiffnessMatrixElement[m, n] += stiffnessMatrixGaussPoint[m, n];
					}
				}
			}
			return Matrix.CreateFromArray(stiffnessMatrixElement.Data);
		}

		private static Matrix2D CalculateDeformationMatrix2(TSplineElement2D tsplineElement,
			ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			Matrix2D B2 = new Matrix2D(4, 2 * tsplineElement.ControlPoints.Count);
			for (int column = 0; column < 2 * tsplineElement.ControlPoints.Count; column += 2)
			{
				B2[0, column] += tsplines.TSplineDerivativeValuesKsi[column / 2, j];
				B2[1, column] += tsplines.TSplineDerivativeValuesHeta[column / 2, j];
				B2[2, column + 1] += tsplines.TSplineDerivativeValuesKsi[column / 2, j];
				B2[3, column + 1] += tsplines.TSplineDerivativeValuesHeta[column / 2, j];
			}

			return B2;
		}

		private static Matrix2D CalculateDeformationMatrix1(Matrix2D jacobianMatrix, double jacdet)
		{
			Matrix2D B1 = new Matrix2D(3, 4);

			B1[0, 0] += jacobianMatrix[1, 1] / jacdet;
			B1[0, 1] += -jacobianMatrix[0, 1] / jacdet;
			B1[1, 2] += -jacobianMatrix[1, 0] / jacdet;
			B1[1, 3] += jacobianMatrix[0, 0] / jacdet;
			B1[2, 0] += -jacobianMatrix[1, 0] / jacdet;
			B1[2, 1] += jacobianMatrix[0, 0] / jacdet;
			B1[2, 2] += jacobianMatrix[1, 1] / jacdet;
			B1[2, 3] += -jacobianMatrix[0, 1] / jacdet;
			return B1;
		}

		private static double CalculateJacobianDeterminant(Matrix2D jacobianMatrix)
		{
			double jacdet = jacobianMatrix[0, 0] * jacobianMatrix[1, 1]
			                - jacobianMatrix[1, 0] * jacobianMatrix[0, 1];
			return jacdet;
		}

		private static Matrix2D CalculateJacobianMatrix(TSplineElement2D tsplineElement,
			ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			Matrix2D jacobianMatrix = new Matrix2D(2, 2);

			for (int k = 0; k < tsplineElement.ControlPoints.Count; k++)
			{
				jacobianMatrix[0, 0] += tsplines.TSplineDerivativeValuesKsi[k, j] * tsplineElement.ControlPoints[k].X;
				jacobianMatrix[0, 1] += tsplines.TSplineDerivativeValuesKsi[k, j] * tsplineElement.ControlPoints[k].Y;
				jacobianMatrix[1, 0] += tsplines.TSplineDerivativeValuesHeta[k, j] * tsplineElement.ControlPoints[k].X;
				jacobianMatrix[1, 1] += tsplines.TSplineDerivativeValuesHeta[k, j] * tsplineElement.ControlPoints[k].Y;
			}

			return jacobianMatrix;
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(TSplineElement2D element)
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
			var tsplineElement = (TSplineElement2D)element;
			var knotParametricCoordinatesKsi = new Vector(new double[] { -1, 1 });
			var knotParametricCoordinatesHeta = new Vector(new double[] { -1, 1 });

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(tsplineElement, tsplineElement.ControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);

			var knotDisplacements = new double[4, 3];
			var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			for (int j = 0; j < knotDisplacements.GetLength(0); j++)
			{
				for (int i = 0; i < element.ControlPoints.Count; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += tsplines.TSplineValues[i, j] * localDisplacements[i, 0];
					knotDisplacements[paraviewKnotRenumbering[j], 1] += tsplines.TSplineValues[i, j] * localDisplacements[i, 1];
					knotDisplacements[paraviewKnotRenumbering[j], 2] += tsplines.TSplineValues[i, j] * localDisplacements[i, 2];
				}
			}

			return knotDisplacements;
		}

		public double[,] CalculatePointsForPostProcessing(TSplineElement2D element)
		{
			var localCoordinates = new double[4, 2]
			{
				{-1, -1},
				{-1, 1},
				{1, -1},
				{1, 1}
			};

			var knotParametricCoordinatesKsi = new Vector(new double[] { -1, 1 });
			var knotParametricCoordinatesHeta = new Vector(new double[] { -1, 1 });

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(element, element.ControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);

			var knotDisplacements = new double[4, 3];
			var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			for (int j = 0; j < localCoordinates.GetLength(0); j++)
			{
				for (int i = 0; i < element.ControlPoints.Count; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += tsplines.TSplineValues[i, j] * element.ControlPoints[i].X;
					knotDisplacements[paraviewKnotRenumbering[j], 1] += tsplines.TSplineValues[i, j] * element.ControlPoints[i].Y;
					knotDisplacements[paraviewKnotRenumbering[j], 2] += tsplines.TSplineValues[i, j] * element.ControlPoints[i].Z;
				}
			}

			return knotDisplacements;

		}
	}
}
