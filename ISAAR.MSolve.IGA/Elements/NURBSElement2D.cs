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
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Problems.Structural.Elements
{
    public class NURBSElement2D : Element, IStructuralIsogeometricElement
	{
		protected readonly static IDofType[] controlPointDOFTypes = new IDofType[] {StructuralDof.TranslationX, StructuralDof.TranslationY};
		protected IDofType[][] dofTypes;
		protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private DynamicMaterial dynamicProperties;
        public CellType CellType { get; } = CellType.Unknown;

        #region IStructuralIsogeometricElement


        public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }

			set { this.dofEnumerator = value; }
		}

		public ElementDimensions ElementDimensions
		{
			get { return ElementDimensions.TwoD; }
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			var nurbsElement = (NURBSElement2D) element;
			dofTypes = new IDofType[nurbsElement.ControlPoints.Count][];
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

		public IMatrix DampingMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IMatrix MassMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IMatrix StiffnessMatrix(IElement element)
		{
			var nurbsElement = (NURBSElement2D) element;
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(nurbsElement);
			Matrix stiffnessMatrixElement = Matrix.CreateZero(nurbsElement.ControlPointsDictionary.Count * 2,
				nurbsElement.ControlPointsDictionary.Count * 2);

			NURBS2D nurbs = new NURBS2D(nurbsElement, nurbsElement.ControlPoints);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				Matrix jacobianMatrix = Matrix.CreateZero(2, 2);

				for (int k = 0; k < nurbsElement.ControlPoints.Count; k++)
				{
					jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * nurbsElement.ControlPoints[k].X;
					jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * nurbsElement.ControlPoints[k].Y;
					jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * nurbsElement.ControlPoints[k].X;
					jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * nurbsElement.ControlPoints[k].Y;
				}

				double jacdet = jacobianMatrix[0, 0] * jacobianMatrix[1, 1]
				                - jacobianMatrix[1, 0] * jacobianMatrix[0, 1];

				Matrix B1 = Matrix.CreateZero(3, 4);

				B1[0, 0] += jacobianMatrix[1, 1] / jacdet;
				B1[0, 1] += -jacobianMatrix[0, 1] / jacdet;
				B1[1, 2] += -jacobianMatrix[1, 0] / jacdet;
				B1[1, 3] += jacobianMatrix[0, 0] / jacdet;
				B1[2, 0] += -jacobianMatrix[1, 0] / jacdet;
				B1[2, 1] += jacobianMatrix[0, 0] / jacdet;
				B1[2, 2] += jacobianMatrix[1, 1] / jacdet;
				B1[2, 3] += -jacobianMatrix[0, 1] / jacdet;

				Matrix B2 = Matrix.CreateZero(4, 2 * nurbsElement.ControlPoints.Count);
				for (int column = 0; column < 2 * nurbsElement.ControlPoints.Count; column += 2)
				{
					B2[0, column] += nurbs.NurbsDerivativeValuesKsi[column / 2, j];
					B2[1, column] += nurbs.NurbsDerivativeValuesHeta[column / 2, j];
					B2[2, column + 1] += nurbs.NurbsDerivativeValuesKsi[column / 2, j];
					B2[3, column + 1] += nurbs.NurbsDerivativeValuesHeta[column / 2, j];
				}

				Matrix B = B1 * B2;
				IMatrixView ElasticityMatrix = ((IContinuumMaterial2D)nurbsElement.Patch.Material).ConstitutiveMatrix;
				Matrix stiffnessMatrixGaussPoint = B.ThisTransposeTimesOtherTimesThis(ElasticityMatrix);
				stiffnessMatrixGaussPoint = stiffnessMatrixGaussPoint *
				                            (jacdet * gaussPoints[j].WeightFactor * nurbsElement.Patch.Thickness);

				for (int m = 0; m < nurbsElement.ControlPoints.Count * 2; m++)
				{
					for (int n = 0; n < nurbsElement.ControlPoints.Count * 2; n++)
					{
						stiffnessMatrixElement[m, n] += stiffnessMatrixGaussPoint[m, n];
					}
				}
			}

			return stiffnessMatrixElement;
		}


		public double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements)
		{
			var nurbsElement = (NURBSElement2D)element;
			var knotParametricCoordinatesKsi= Vector.CreateFromArray(new double[]{element.Knots[0].Ksi, element.Knots[2].Ksi });
			var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] { element.Knots[0].Heta, element.Knots[1].Heta });
			NURBS2D nurbs = new NURBS2D(nurbsElement, nurbsElement.ControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);
			var knotDisplacements = new double[4, 2];
			var paraviewKnotRenumbering = new int[] {0, 3, 1, 2};
			for (int j = 0; j < element.Knots.Count; j++)
			{
				for (int i = 0; i < element.ControlPoints.Count; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
					knotDisplacements[paraviewKnotRenumbering[j], 1] += nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
				}
			}

			return knotDisplacements;
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
			NeumannBoundaryCondition neumann)
		{
			throw new NotSupportedException();
		}


		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
			NeumannBoundaryCondition neumann)
		{
			IList<GaussLegendrePoint3D> gaussPoints =
				CreateElementGaussPoints(element, face.Degrees[0], face.Degrees[1]);
			Dictionary<int, double> neumannLoad = new Dictionary<int, double>();

			NURBS2D nurbs = new NURBS2D(element, element.ControlPoints, face);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
				double xGaussPoint = 0;
				double yGaussPoint = 0;
				double zGaussPoint = 0;
				for (int k = 0; k < element.ControlPoints.Count; k++)
				{
					xGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].X;
					yGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].Y;
					zGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].Z;
					jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].X;
					jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].Y;
					jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].Z;
					jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * element.ControlPoints[k].X;
					jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * element.ControlPoints[k].Y;
					jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * element.ControlPoints[k].Z;
				}

				Vector surfaceBasisVector1 = Vector.CreateZero(3);
				surfaceBasisVector1[0] = jacobianMatrix[0, 0];
				surfaceBasisVector1[1] = jacobianMatrix[0, 1];
				surfaceBasisVector1[2] = jacobianMatrix[0, 2];

				Vector surfaceBasisVector2 = Vector.CreateZero(3);
				surfaceBasisVector2[0] = jacobianMatrix[1, 0];
				surfaceBasisVector2[1] = jacobianMatrix[1, 1];
				surfaceBasisVector2[2] = jacobianMatrix[1, 2];

				Vector surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);

				double jacdet = jacobianMatrix[0, 0] * jacobianMatrix[1, 1]
				                - jacobianMatrix[1, 0] * jacobianMatrix[0, 1];

				for (int k = 0; k < element.ControlPoints.Count; k++)
				{
					//int dofIDX = element.Model.ControlPointDOFsDictionary[element.ControlPoints[k].ID][DOFType.X];
					//int dofIDY = element.Model.ControlPointDOFsDictionary[element.ControlPoints[k].ID][DOFType.Y];
					//int dofIDZ = element.Model.ControlPointDOFsDictionary[element.ControlPoints[k].ID][DOFType.Y];
					int dofIDX = element.Model.GlobalDofOrdering.GlobalFreeDofs[element.ControlPoints[k], StructuralDof.TranslationX];
					int dofIDY = element.Model.GlobalDofOrdering.GlobalFreeDofs[element.ControlPoints[k], StructuralDof.TranslationY];
					int dofIDZ = element.Model.GlobalDofOrdering.GlobalFreeDofs[element.ControlPoints[k], StructuralDof.TranslationZ];

					if (neumannLoad.ContainsKey(dofIDX))
						neumannLoad[dofIDX] += nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
						                       neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[0] *
						                       surfaceBasisVector3[0];
					else
						neumannLoad.Add(dofIDX,
							nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
							neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[0] * surfaceBasisVector3[0]);

					if (neumannLoad.ContainsKey(dofIDY))
						neumannLoad[dofIDY] += nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
						                       neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[1] *
						                       surfaceBasisVector3[1];
					else
						neumannLoad.Add(dofIDY,
							nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
							neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[1] * surfaceBasisVector3[1]);

					if (neumannLoad.ContainsKey(dofIDZ))
						neumannLoad[dofIDZ] += nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
						                       neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[2] *
						                       surfaceBasisVector3[2];
					else
						neumannLoad.Add(dofIDZ,
							nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
							neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[2] * surfaceBasisVector3[2]);
				}
			}

			return neumannLoad;
		}

		#endregion

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
			PressureBoundaryCondition pressure)
		{
			var dofs = new IDofType[] {StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ};

			IList<GaussLegendrePoint3D> gaussPoints =
				CreateElementGaussPoints(element, face.Degrees[0], face.Degrees[1]);
			Dictionary<int, double> pressureLoad = new Dictionary<int, double>();

			NURBS2D nurbs = new NURBS2D(element, element.ControlPoints, face);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
				double xGaussPoint = 0;
				double yGaussPoint = 0;
				double zGaussPoint = 0;
				for (int k = 0; k < element.ControlPoints.Count; k++)
				{
					xGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].X;
					yGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].Y;
					zGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].Z;
					jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].X;
					jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].Y;
					jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].Z;
					jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * element.ControlPoints[k].X;
					jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * element.ControlPoints[k].Y;
					jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * element.ControlPoints[k].Z;
				}

				Vector surfaceBasisVector1 = Vector.CreateZero(3);
				surfaceBasisVector1[0] = jacobianMatrix[0, 0];
				surfaceBasisVector1[1] = jacobianMatrix[0, 1];
				surfaceBasisVector1[2] = jacobianMatrix[0, 2];

				Vector surfaceBasisVector2 = Vector.CreateZero(3);
				surfaceBasisVector2[0] = jacobianMatrix[1, 0];
				surfaceBasisVector2[1] = jacobianMatrix[1, 1];
				surfaceBasisVector2[2] = jacobianMatrix[1, 2];

				Vector surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);

				double jacdet = jacobianMatrix[0, 0] * jacobianMatrix[1, 1]
				                - jacobianMatrix[1, 0] * jacobianMatrix[0, 1];

				for (int k = 0; k < element.ControlPoints.Count; k++)
				{
					for (int m = 0; m < 3; m++)
					{
						int dofID = element.Model.GlobalDofOrdering.GlobalFreeDofs[element.ControlPoints[k],dofs[m]];
						;
						if (pressureLoad.ContainsKey(dofID))
						{
							pressureLoad[dofID] += nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
							                       pressure.Value * surfaceBasisVector3[m];
						}
						else
						{
							pressureLoad.Add(dofID,
								nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor * pressure.Value *
								surfaceBasisVector3[m]);
						}
					}
				}
			}

			return pressureLoad;
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
			PressureBoundaryCondition pressure)
		{
			throw new NotSupportedException();
		}


		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Knots);
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element, int degreeKsi, int degreeHeta)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, element.Knots);
		}

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}
	}
}