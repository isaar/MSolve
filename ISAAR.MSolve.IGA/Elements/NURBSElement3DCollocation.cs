using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
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
using ISAAR.MSolve.Solvers.LinearSystems;
using Vector = ISAAR.MSolve.LinearAlgebra.Vectors.Vector;

namespace ISAAR.MSolve.IGA.Elements
{
	/// <summary>
	/// Three dimensional collocation point.
	/// Calculated according to "Improving the Computational Performance of Isogeometric Analysis" of Gkritzalis, Christos.
	/// Authors: Dimitris Tsapetis.
	/// </summary>
	public class NURBSElement3DCollocation:Element, IStructuralIsogeometricElement, ICollocationElement
	{
        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
        protected IDofType[][] dofTypes;
        protected readonly static IDofType[] controlPointDOFTypes = new StructuralDof[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
        public CellType CellType { get; } = CellType.Unknown;


        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;
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
        public bool MaterialModified { get; }

        private CollocationPoint3D _collocationPoint;
        public CollocationPatch Patch { get; set; }

        public new IStructuralAsymmetricModel Model { get; set; }
        public CollocationPoint3D CollocationPoint
        {
            get => _collocationPoint;
            set => _collocationPoint = value;
        }

        public IList<IDofType> GetDOFTypesForDOFEnumeration(IElement element)
        {
            return new StructuralDof[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
        }

        INode ICollocationElement.CollocationPoint { get => _collocationPoint; set => _collocationPoint = (CollocationPoint3D)value; }
        IAsymmetricSubdomain ICollocationElement.Patch { get => this.Patch; set => this.Patch=(CollocationPatch)value; }
       

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

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
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
			var elementCollocation = (NURBSElement3DCollocation)element;

            var nurbs = new NURBS3D(elementCollocation.Patch.NumberOfControlPointsKsi,
                elementCollocation.Patch.NumberOfControlPointsHeta, elementCollocation.Patch.NumberOfControlPointsZeta,
                elementCollocation.Patch.DegreeKsi, elementCollocation.Patch.DegreeHeta,
                elementCollocation.Patch.DegreeZeta, elementCollocation.Patch.KnotValueVectorKsi,
                elementCollocation.Patch.KnotValueVectorHeta, elementCollocation.Patch.KnotValueVectorZeta,
                elementCollocation.ControlPoints.ToArray(), elementCollocation.CollocationPoint);

            var jacobianMatrix = CalculateJacobian(elementCollocation.ControlPoints.ToArray(), nurbs, 0);

            var inverseJacobian = jacobianMatrix.Invert();

            var dR = CalculateNaturalDerivatives(nurbs, inverseJacobian);

            if (elementCollocation.CollocationPoint.IsBoundary)
            {
                var (normalX, normalY, normalZ) = CalculateNormalVectors(elementCollocation, jacobianMatrix);
                return CalculateCollocationPointStiffnessBoundary(elementCollocation, normalX, normalY, normalZ, dR);
            }
            else
            {
                var hessianMatrix = CalculateHessian(elementCollocation.ControlPoints, nurbs,0);
                var squareDerivatives = CalculateSquareDerivatives(jacobianMatrix);
                var ddR = CalculateNaturalSecondDerivatives(nurbs, hessianMatrix, dR, squareDerivatives);

                return CollocationPointStiffness(elementCollocation, ddR);
            }
			
		}


        private Matrix CalculateCollocationPointStiffnessBoundary(NURBSElement3DCollocation elementCollocation,
            double xGaussPoint, double yGaussPoint, double zGaussPoint, Matrix dR )
        {
            var collocationPointStiffness = Matrix.CreateZero(3, elementCollocation.ControlPoints.Count * 3);
            var E = elementCollocation.Patch.Material.YoungModulus;
            var nu = elementCollocation.Patch.Material.PoissonRatio;
            var lambda = E * nu / (1 + nu) / (1 - 2 * nu);
            var m = E / (2 * (1 + nu));

            for (int i = 0; i < elementCollocation.ControlPoints.Count * 3; i += 3)
            {
                var index = i / 3;
                collocationPointStiffness[0, i] = (lambda + 2 * m) * xGaussPoint * dR[0, index] +
                                                  m * yGaussPoint * dR[1, index] + m * zGaussPoint * dR[2, index];
                collocationPointStiffness[0, i + 1] = lambda * xGaussPoint * dR[1, index] + m * yGaussPoint * dR[0, index];
                collocationPointStiffness[0, i + 2] = lambda * xGaussPoint * dR[2, index] + m * zGaussPoint * dR[0, index];

                collocationPointStiffness[1, i] = lambda * yGaussPoint * dR[0, index] + m * xGaussPoint * dR[1, index];
                collocationPointStiffness[1, i + 1] = (lambda + 2 * m) * yGaussPoint * dR[1, index] +
                                                      m * xGaussPoint * dR[0, index] + m * zGaussPoint * dR[2, index];
                collocationPointStiffness[1, i + 2] = lambda * yGaussPoint * dR[2, index] + m * zGaussPoint * dR[1, index];

                collocationPointStiffness[2, i] = lambda * zGaussPoint * dR[0, index] + m * xGaussPoint * dR[2, index];
                collocationPointStiffness[2, i + 1] = lambda * zGaussPoint * dR[1, index] + m * yGaussPoint * dR[2, index];
                collocationPointStiffness[2, i + 2] = (lambda + 2 * m) * zGaussPoint * dR[2, index] +
                                                      m * xGaussPoint * dR[0, index] + m * yGaussPoint * dR[1, index];
            }

            return collocationPointStiffness;
        }

        private (double xGaussPoint, double yGaussPoint, double zGaussPoint) CalculateNormalVectors(
            NURBSElement3DCollocation elementCollocation, Matrix3by3 jacobianMatrix)
        {
            List<Vector> normals= new List<Vector>();
            foreach (var surface in elementCollocation.CollocationPoint.Surfaces)
            {
                var u = Vector.CreateFromArray(new double[] { jacobianMatrix[0, 0], jacobianMatrix[0, 1], jacobianMatrix[0, 2] });
                var v = Vector.CreateFromArray(new double[] { jacobianMatrix[1, 0], jacobianMatrix[1, 1], jacobianMatrix[1, 2] });
                var w = Vector.CreateFromArray(new double[] { jacobianMatrix[2, 0], jacobianMatrix[2, 1], jacobianMatrix[2, 2] });

                if (surface == Surface.Bottom)
                {
                    var normal = v.CrossProduct(u);
                    normals.Add(normal.Scale(1/normal.Norm2()));
                }
                else if (surface == Surface.Top)
                {
                    var normal = u.CrossProduct(v);
                    normals.Add(normal.Scale(1 / normal.Norm2()));
                }
                else if (surface == Surface.Front)
                {
                    var normal = u.CrossProduct(w);
                    normals.Add(normal.Scale(1 / normal.Norm2()));
                }
                else if (surface == Surface.Back)
                {
                    var normal = w.CrossProduct(u);
                    normals.Add(normal.Scale(1 / normal.Norm2()));
                }
                else if (surface==Surface.Left)
                {
                    var normal = w.CrossProduct(v);
                    normals.Add(normal.Scale(1 / normal.Norm2()));
                }
                else if (surface == Surface.Right)
                {
                    var normal = w.CrossProduct(v);
                    normals.Add(normal.Scale(1 / normal.Norm2()));
                }
            }

            var normalVector = Vector3.CreateZero();
            normalVector[0] = normals.Sum(n => n[0]);
            normalVector[1] = normals.Sum(n => n[1]);
            normalVector[2] = normals.Sum(n => n[2]);
            return (normalVector[0], normalVector[1], normalVector[2]);
        }

        private static Matrix CollocationPointStiffness(NURBSElement3DCollocation elementCollocation, Matrix ddR)
		{
			var collocationPointStiffness = Matrix.CreateZero(3, elementCollocation.ControlPoints.Count * 3);

			var E = elementCollocation.Patch.Material.YoungModulus;
			var nu = elementCollocation.Patch.Material.PoissonRatio;

			var lambda = nu * E / (1 + nu) / (1 - 2 * nu);
			var m = E / 2 / (1 + nu);

			for (int i = 0; i < elementCollocation.ControlPoints.Count * 3; i += 3)
			{
				var index = i / 3;
				collocationPointStiffness[0, i] = (lambda + 2 * m) * ddR[0, index] + m * ddR[1, index] + m * ddR[2, index];
				collocationPointStiffness[1, i] = (lambda + m) * ddR[3, index];
				collocationPointStiffness[2, i] = (lambda + m) * ddR[4, index];

				collocationPointStiffness[0, i + 1] = (lambda + m) * ddR[3, index];
				collocationPointStiffness[1, i + 1] = (lambda + 2 * m) * ddR[1, index] + m * ddR[0, index] + m * ddR[2, index];
				collocationPointStiffness[2, i + 1] = (lambda + m) * ddR[5, index];

				collocationPointStiffness[0, i + 2] = (lambda + m) * ddR[4, index];
				collocationPointStiffness[1, i + 2] = (lambda + m) * ddR[5, index];
				collocationPointStiffness[2, i + 2] = (lambda + 2 * m) * ddR[2, index] + m * ddR[0, index] + m * ddR[1, index];
			}

			return collocationPointStiffness;
		}

		private Matrix CalculateNaturalSecondDerivatives(NURBS3D nurbs, Matrix hessianMatrix, Matrix dR, Matrix squareDerivatives)
		{
			var ddR2 = hessianMatrix * dR;

            var ddR3 = Matrix.CreateZero(6, nurbs.NurbsSecondDerivativeValueKsi.NumRows);
			for (int i = 0; i < ddR3.NumColumns; i++)
			{
				ddR3[0, i] = nurbs.NurbsSecondDerivativeValueKsi[i, 0] - ddR2[0, i];
				ddR3[1, i] = nurbs.NurbsSecondDerivativeValueHeta[i, 0] - ddR2[1, i];
				ddR3[2, i] = nurbs.NurbsSecondDerivativeValueZeta[i, 0] - ddR2[2, i];

				ddR3[3, i] = nurbs.NurbsSecondDerivativeValueKsiHeta[i, 0] - ddR2[3, i];
				ddR3[4, i] = nurbs.NurbsSecondDerivativeValueKsiZeta[i, 0] - ddR2[4, i];
				ddR3[5, i] = nurbs.NurbsSecondDerivativeValueHetaZeta[i, 0] - ddR2[5, i];
			}
            
            var squareInvert = squareDerivatives.Invert();

            var ddR = squareInvert * ddR3;


   //         for (int i = 0; i < ddR.GetLength(1); i++)
			//{
			//	var rhs = Vector.CreateFromArray(new[]
			//		{ddR3[0, i], ddR3[1, i], ddR3[2, i], ddR3[3, i], ddR3[4, i], ddR3[5, i]});
   //             var solution = squareInvert * rhs;
			//	ddR[0, i] = solution[0];
			//	ddR[1, i] = solution[1];
			//	ddR[2, i] = solution[2];
			//	ddR[3, i] = solution[3];
			//	ddR[4, i] = solution[4];
			//	ddR[5, i] = solution[5];
			//}

			return ddR;
		}

		private Matrix CalculateNaturalDerivatives(NURBS3D nurbs, Matrix3by3 inverseJacobian)
		{
			var dR = Matrix.CreateZero(3, nurbs.NurbsDerivativeValuesKsi.NumRows);
			for (int i = 0; i < dR.NumColumns; i++)
			{
				var dKsi = nurbs.NurbsDerivativeValuesKsi[i, 0];
				var dHeta = nurbs.NurbsDerivativeValuesHeta[i, 0];
				var dZeta = nurbs.NurbsDerivativeValuesZeta[i, 0];

				dR[0, i] = inverseJacobian[0, 0] * dKsi + inverseJacobian[0, 1] * dHeta + inverseJacobian[0, 2] * dZeta;
				dR[1, i] = inverseJacobian[1, 0] * dKsi + inverseJacobian[1, 1] * dHeta + inverseJacobian[1, 2] * dZeta;
				dR[2, i] = inverseJacobian[2, 0] * dKsi + inverseJacobian[2, 1] * dHeta + inverseJacobian[2, 2] * dZeta;
			}

			return dR;
		}

		private Matrix CalculateSquareDerivatives(Matrix3by3 jacobianMatrix)
		{
			var J11 = jacobianMatrix[0, 0];
			var J12 = jacobianMatrix[0, 1];
			var J13 = jacobianMatrix[0, 2];
			var J21 = jacobianMatrix[1, 0];
			var J22 = jacobianMatrix[1, 1];
			var J23 = jacobianMatrix[1, 2];
			var J31 = jacobianMatrix[2, 0];
			var J32 = jacobianMatrix[2, 1];
			var J33 = jacobianMatrix[2, 2];


			var squareDerivatives = Matrix.CreateZero(6, 6);
			squareDerivatives[0, 0] = J11*J11;
			squareDerivatives[0, 1] = J12*J12;
			squareDerivatives[0, 2] = J13*J13;
			squareDerivatives[0, 3] = 2*J11*J12;
			squareDerivatives[0, 4] = 2*J11*J13;
			squareDerivatives[0, 5] = 2*J12*J13;

			squareDerivatives[1, 0] = J21*J21;
			squareDerivatives[1, 1] = J22*J22;
			squareDerivatives[1, 2] = J23*J23;
			squareDerivatives[1, 3] = 2*J21*J22;
			squareDerivatives[1, 4] = 2*J21*J23;
			squareDerivatives[1, 5] = 2*J22*J23;

			squareDerivatives[2, 0] = J31*J31;
			squareDerivatives[2, 1] = J32*J32;
			squareDerivatives[2, 2] = J33*J33;
			squareDerivatives[2, 3] = 2*J31*J32;
			squareDerivatives[2, 4] = 2*J31*J33;
			squareDerivatives[2, 5] = 2*J32*J33;

			squareDerivatives[3, 0] = J11*J21;
			squareDerivatives[3, 1] = J12*J22;
			squareDerivatives[3, 2] = J13*J23;
			squareDerivatives[3, 3] = J11*J22+J21*J12;
			squareDerivatives[3, 4] = J11*J23+J21*J13;
			squareDerivatives[3, 5] = J12*J23+J22*J13;

			squareDerivatives[4, 0] = J11*J31;
			squareDerivatives[4, 1] = J12*J32;
			squareDerivatives[4, 2] = J13*J33;
			squareDerivatives[4, 3] = J11*J32+J31*J12;
			squareDerivatives[4, 4] = J11*J33+J31*J13;
			squareDerivatives[4, 5] = J12*J33+J32*J13;

			squareDerivatives[5, 0] = J21*J31;
			squareDerivatives[5, 1] = J22*J32;
			squareDerivatives[5, 2] = J23*J33;
			squareDerivatives[5, 3] = J21*J32+J31*J22;
			squareDerivatives[5, 4] = J21*J33+J31*J23;
			squareDerivatives[5, 5] = J22*J33+J32*J23;

			return squareDerivatives;
		}

		private Matrix CalculateHessian(IList<ControlPoint> controlPoints, NURBS3D nurbs, int j)
		{
			var hessianMatrix = Matrix.CreateZero(6,3);
			for (int k = 0; k < controlPoints.Count; k++)
			{
				hessianMatrix[0, 0] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * controlPoints[k].X;
				hessianMatrix[0, 1] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * controlPoints[k].Y;
				hessianMatrix[0, 2] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * controlPoints[k].Z;

				hessianMatrix[1, 0] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * controlPoints[k].X;
				hessianMatrix[1, 1] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * controlPoints[k].Y;
				hessianMatrix[1, 2] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * controlPoints[k].Z;

				hessianMatrix[2, 0] += nurbs.NurbsSecondDerivativeValueZeta[k, j] * controlPoints[k].X;
				hessianMatrix[2, 1] += nurbs.NurbsSecondDerivativeValueZeta[k, j] * controlPoints[k].Y;
				hessianMatrix[2, 2] += nurbs.NurbsSecondDerivativeValueZeta[k, j] * controlPoints[k].Z;

				hessianMatrix[3, 0] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * controlPoints[k].X;
				hessianMatrix[3, 1] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * controlPoints[k].Y;
				hessianMatrix[3, 2] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * controlPoints[k].Z;

				hessianMatrix[4, 0] += nurbs.NurbsSecondDerivativeValueKsiZeta[k, j] * controlPoints[k].X;
				hessianMatrix[4, 1] += nurbs.NurbsSecondDerivativeValueKsiZeta[k, j] * controlPoints[k].Y;
				hessianMatrix[4, 2] += nurbs.NurbsSecondDerivativeValueKsiZeta[k, j] * controlPoints[k].Z;

				hessianMatrix[5, 0] += nurbs.NurbsSecondDerivativeValueHetaZeta[k, j] * controlPoints[k].X;
				hessianMatrix[5, 1] += nurbs.NurbsSecondDerivativeValueHetaZeta[k, j] * controlPoints[k].Y;
				hessianMatrix[5, 2] += nurbs.NurbsSecondDerivativeValueHetaZeta[k, j] * controlPoints[k].Z;
			}

			return hessianMatrix;
		}

		private Matrix3by3 CalculateJacobian(ControlPoint[] elementControlPoints, NURBS3D nurbs, int j)
		{
			var jacobianMatrix = Matrix3by3.CreateZero();

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
            var numberOfCP= ((NURBSElement3DCollocation)element).ControlPointsDictionary.Count;
            dofTypes = new IDofType[numberOfCP][];
            for (int i = 0; i < numberOfCP; i++)
            {
                dofTypes[i] = controlPointDOFTypes;
            }

            return dofTypes;
        }
	}
}
