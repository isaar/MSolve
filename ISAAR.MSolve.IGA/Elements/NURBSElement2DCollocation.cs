using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.IGA.Elements
{
	public class NURBSElement2DCollocation : Element, IStructuralIsogeometricElement, ICollocationElement
    {
        protected readonly static IDofType[] controlPointDOFTypes = new StructuralDof[] { StructuralDof.TranslationX, StructuralDof.TranslationY };
        protected IDofType[][] dofTypes;
        private CollocationPoint2D _collocationPoint;
        public CellType CellType { get; } = CellType.Unknown;

        public CollocationPoint2D CollocationPoint
        {
            get => _collocationPoint;
            set => _collocationPoint = value;
        }
        public CollocationPatch Patch { get; set; }
        public IStructuralAsymmetricModel Model { get; set; }

        public IList<IDofType> GetDOFTypesForDOFEnumeration(IElement element)
        {
            return new StructuralDof[] { StructuralDof.TranslationX, StructuralDof.TranslationY };
        }

        IAsymmetricSubdomain ICollocationElement.Patch
        {
            get => Patch;
            set => Patch = (CollocationPatch)value;
        }

        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;
        public IElementDofEnumerator DofEnumerator
        {
            get { return dofEnumerator; }

            set { this.dofEnumerator = value; }
        }
        public bool MaterialModified { get; }
        INode ICollocationElement.CollocationPoint { get => _collocationPoint; set => _collocationPoint=(CollocationPoint2D)value; }

        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
			NeumannBoundaryCondition neumann)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
			NeumannBoundaryCondition neumann)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
			PressureBoundaryCondition pressure)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
			PressureBoundaryCondition pressure)
		{
			throw new NotImplementedException();
		}

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements,
			double[] localdDisplacements)
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

		public double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public void ClearMaterialState()
		{
			throw new NotImplementedException();
		}

		public IMatrix StiffnessMatrix(IElement element)
		{

			var elementCollocation = (NURBSElement2DCollocation) element;

            var nurbs = new NURBS2D(elementCollocation.Patch.DegreeKsi, elementCollocation.Patch.DegreeHeta,
                elementCollocation.Patch.KnotValueVectorKsi, elementCollocation.Patch.KnotValueVectorHeta,
                elementCollocation.CollocationPoint, elementCollocation.ControlPoints);

            var jacobianMatrix = CalculateJacobianMatrix(elementCollocation, nurbs);

            var inverseJacobian = jacobianMatrix.Invert();
            var dR = CalculateNaturalDerivatives(nurbs, inverseJacobian);
            //if (elementCollocation.Patch.Material is ElasticMaterial2D elasticMaterial2D&& elasticMaterial2D.StressState!=StressState2D.PlaneStress)
            //    throw  new NotSupportedException("Cannot use Plane stress with collocation for the time being");

            if (elementCollocation.CollocationPoint.IsBoundary)
            {
                var (xGaussPoint, yGaussPoint) = CalculateNormalVectors(elementCollocation, nurbs);

                return CalculateCollocationPointStiffnessBoundary(elementCollocation, xGaussPoint, dR, yGaussPoint);
            }
            else
            {
                var hessianMatrix = CalculateHessian(elementCollocation, nurbs, 0);
                var squareDerivatives = CalculateSquareDerivatives(jacobianMatrix);
                var ddR = CalculateNaturalSecondDerivatives(nurbs, hessianMatrix, dR, squareDerivatives);

                return CalculateCollocationPointStiffness(elementCollocation, ddR);
            }
			
		}

        private Matrix CalculateCollocationPointStiffnessBoundary(NURBSElement2DCollocation elementCollocation,
            double xGaussPoint, double[,] dR, double yGaussPoint)
        {
            var collocationPointStiffness = Matrix.CreateZero(2, elementCollocation.ControlPoints.Count * 2);
            var E = elementCollocation.Patch.Material.YoungModulus;
            var nu = elementCollocation.Patch.Material.PoissonRatio;
            //var lambda = E * nu / (1 + nu) / (1 - 2 * nu);
            //var m = E / (2 * (1 + nu));
            var aux = E / (1 - nu * nu);

            for (int i = 0; i < elementCollocation.ControlPoints.Count * 2; i += 2)
            {
                var index = i / 2;
                collocationPointStiffness[0, i] = aux * (xGaussPoint * dR[0, index] + yGaussPoint * (1 - nu) / 2 * dR[1, index]);
                collocationPointStiffness[0, i + 1] = aux*(xGaussPoint*nu*dR[1,index]+yGaussPoint*(1-nu)/2*dR[0,index]);

                collocationPointStiffness[1, i] = aux * (yGaussPoint * nu * dR[0, index] + xGaussPoint * (1 - nu) / 2 * dR[1, index]);
                collocationPointStiffness[1, i + 1] = aux * (yGaussPoint * dR[1, index] + xGaussPoint * (1 - nu) / 2 * dR[0, index]);
            }

            return collocationPointStiffness;
        }

        private (double xGaussPoint, double yGaussPoint) CalculateNormalVectors(
            NURBSElement2DCollocation elementCollocation, NURBS2D nurbs)
        {
            double xGaussPoint = 0;
            double yGaussPoint = 0;
            for (int k = 0; k < elementCollocation.ControlPoints.Count; k++)
            {
                xGaussPoint += nurbs.NurbsValues[k, 0] * elementCollocation.ControlPoints[k].X;
                yGaussPoint += nurbs.NurbsValues[k, 0] * elementCollocation.ControlPoints[k].Y;
            }

            double norm = Math.Sqrt(Math.Pow(xGaussPoint, 2) + Math.Pow(yGaussPoint, 2));
            xGaussPoint /= norm;
            yGaussPoint /= norm;
            return (xGaussPoint, yGaussPoint);
        }


        public Matrix3by3 CalculateSquareDerivatives(Matrix2by2 jacobianMatrix)
		{
			var squareDerivatives = Matrix3by3.CreateFromArray(new double[3, 3]
			{
				{
					jacobianMatrix[0, 0] * jacobianMatrix[0, 0], jacobianMatrix[0, 1] * jacobianMatrix[0, 1], 2 * jacobianMatrix[0, 0] * jacobianMatrix[0, 1]
				},
				{
					jacobianMatrix[1, 0] * jacobianMatrix[1, 0],
					jacobianMatrix[1, 1] * jacobianMatrix[1, 1],
					2 * jacobianMatrix[1, 0] * jacobianMatrix[1, 1]
				},
				{
					jacobianMatrix[0, 0] * jacobianMatrix[1, 0],
					jacobianMatrix[0, 1] * jacobianMatrix[1, 1],
					jacobianMatrix[0, 0] * jacobianMatrix[1, 1] + jacobianMatrix[0, 1] * jacobianMatrix[1, 0],
				},
			}, false);
			return squareDerivatives;
		}

		public Matrix CalculateCollocationPointStiffness(NURBSElement2DCollocation elementCollocation, double[,] ddR)
		{
			var collocationPointStiffness = Matrix.CreateZero(2, elementCollocation.ControlPoints.Count * 2);

			var E = elementCollocation.Patch.Material.YoungModulus;
			var nu = elementCollocation.Patch.Material.PoissonRatio;
			var temp = E / (1 - nu*nu);
			for (int i = 0; i < elementCollocation.ControlPoints.Count * 2; i += 2)
			{
				var index = i / 2;
				collocationPointStiffness[0, i] = (ddR[0, index] + (1 - nu) / 2 * ddR[1, index]) * temp;
				collocationPointStiffness[1, i] = (1 + nu) / 2 * ddR[2, index] * temp;

				collocationPointStiffness[0, i + 1] = (1 + nu) / 2 * ddR[2, index] * temp;
				collocationPointStiffness[1, i + 1] = (ddR[1, index] + (1 - nu) / 2 * ddR[0, index]) * temp;
			}

			return collocationPointStiffness;
		}

		public double[,] CalculateNaturalSecondDerivatives(NURBS2D nurbs, Matrix hessianMatrix, double[,] dR,
			Matrix3by3 squareDerivatives)
		{
			var ddR2 = new double[3, nurbs.NurbsSecondDerivativeValueKsi.NumRows];
			for (int i = 0; i < ddR2.GetLength(1); i++)
			{
				ddR2[0, i] = hessianMatrix[0, 0] * dR[0, i] + hessianMatrix[0, 1] * dR[1, i];
				ddR2[1, i] = hessianMatrix[1, 0] * dR[0, i] + hessianMatrix[1, 1] * dR[1, i];
				ddR2[2, i] = hessianMatrix[2, 0] * dR[0, i] + hessianMatrix[2, 1] * dR[1, i];
			}

			var ddR = new double[3, nurbs.NurbsSecondDerivativeValueKsi.NumRows];
			var squareInvert = squareDerivatives.Invert();
			
			var ddR3 = new double[3, nurbs.NurbsSecondDerivativeValueKsi.NumRows];
			for (int i = 0; i < ddR2.GetLength(1); i++)
			{
				ddR3[0, i] = nurbs.NurbsSecondDerivativeValueKsi[i, 0] - ddR2[0, i];
				ddR3[1, i] = nurbs.NurbsSecondDerivativeValueHeta[i, 0] - ddR2[1, i];
				ddR3[2, i] = nurbs.NurbsSecondDerivativeValueKsiHeta[i, 0] - ddR2[2, i];
			}

			for (int i = 0; i < ddR.GetLength(1); i++)
			{
				ddR[0, i] = squareInvert[0, 0] * ddR3[0, i] + squareInvert[0, 1] * ddR3[1, i] + squareInvert[0, 2] * ddR3[2, i];
				ddR[1, i] = squareInvert[1, 0] * ddR3[0, i] + squareInvert[1, 1] * ddR3[1, i] + squareInvert[1, 2] * ddR3[2, i];
				ddR[2, i] = squareInvert[2, 0] * ddR3[0, i] + squareInvert[2, 1] * ddR3[1, i] + squareInvert[2, 2] * ddR3[2, i];
			}

			return ddR;
		}

		public double[,] CalculateNaturalDerivatives(NURBS2D nurbs, Matrix2by2 inverseJacobian)
		{
			var dR = new double[2, nurbs.NurbsSecondDerivativeValueKsi.NumRows];
			for (int i = 0; i < dR.GetLength(1); i++)
			{
				var dKsi = nurbs.NurbsDerivativeValuesKsi[i, 0];
				var dHeta = nurbs.NurbsDerivativeValuesHeta[i, 0];

				dR[0, i] = inverseJacobian[0, 0] * dKsi + inverseJacobian[0, 1] * dHeta;
				dR[1, i] = inverseJacobian[1, 0] * dKsi + inverseJacobian[1, 1] * dHeta;
			}

			return dR;
		}

		public Matrix2by2 CalculateJacobianMatrix(NURBSElement2DCollocation elementCollocation, NURBS2D nurbs)
		{
			var jacobianMatrix = Matrix2by2.CreateZero();
			for (int k = 0; k < elementCollocation.ControlPoints.Count; k++)
			{
				jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, 0] * elementCollocation.ControlPoints[k].X;
				jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, 0] * elementCollocation.ControlPoints[k].Y;
				jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, 0] * elementCollocation.ControlPoints[k].X;
				jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, 0] * elementCollocation.ControlPoints[k].Y;
			}

			return jacobianMatrix;
		}

		public Vector2 CalculateCartesianCollocationPoint(NURBSElement2DCollocation elementCollocation, NURBS2D nurbs)
		{
			var cartesianCollocationPoint = Vector2.CreateZero();
			for (int k = 0; k < elementCollocation.ControlPoints.Count; k++)
			{
				cartesianCollocationPoint[0] += nurbs.NurbsValues[k, 0] * elementCollocation.ControlPoints[k].X;
				cartesianCollocationPoint[1] += nurbs.NurbsValues[k, 0] * elementCollocation.ControlPoints[k].Y;
			}

			return cartesianCollocationPoint;
		}

		public Matrix CalculateHessian(NURBSElement2DCollocation shellElement, NURBS2D nurbs, int j)
		{
			Matrix hessianMatrix = Matrix.CreateZero(3,2);
			for (int k = 0; k < shellElement.ControlPoints.Count; k++)
			{
				hessianMatrix[0, 0] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[0, 1] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[1, 0] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[1, 1] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[2, 0] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[2, 1] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * shellElement.ControlPoints[k].Y;
			}

			return hessianMatrix;
		}

		public IMatrix MassMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IMatrix DampingMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

        public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
        {
            var nurbsElement = (NURBSElement2DCollocation)element;
            dofTypes = new IDofType[nurbsElement.ControlPoints.Count][];
            for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
            {
                dofTypes[i] = controlPointDOFTypes;
            }

            return dofTypes;
        }
	}
}