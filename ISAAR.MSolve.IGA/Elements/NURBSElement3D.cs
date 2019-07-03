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
    public class NURBSElement3D : Element, IStructuralIsogeometricElement
    {
        protected readonly static IDofType[] controlPointDOFTypes = new  IDofType[] {StructuralDof.TranslationX, StructuralDof.TranslationY , StructuralDof.TranslationZ };
        protected IDofType[][] dofTypes;
        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
	    private DynamicMaterial dynamicProperties;
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

        public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
        {
	        var nurbsElement = (NURBSElement3D) element;
			dofTypes = new IDofType[nurbsElement.ControlPoints.Count][];
            for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
            {
                dofTypes[i] = controlPointDOFTypes;
            }
            return dofTypes;
        }

        public bool MaterialModified
        {
            get
            {
                throw new NotImplementedException();
            }
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

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
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
	        var nurbsElement = (NURBSElement3D)element;
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(nurbsElement);
			Matrix stiffnessMatrixElement = Matrix.CreateZero(nurbsElement.ControlPointsDictionary.Count * 3, nurbsElement.ControlPointsDictionary.Count * 3);

            NURBS3D nurbs = new NURBS3D(nurbsElement, nurbsElement.ControlPoints);

            for (int j = 0; j < gaussPoints.Count; j++)
            {
                Matrix jacobianMatrix = CalculateJacobian(nurbsElement.ControlPoints, nurbs, j);

                double jacdet = CalculateJacobianDeterminant(jacobianMatrix);

                Matrix inverseJacobian = CalculateInverseJacobian(jacobianMatrix, jacdet);

                Matrix B1 = CalculateDeformationMatrix1(inverseJacobian);

                Matrix B2 = CalculateDeformationMatrix2(nurbsElement.ControlPoints, nurbs, j);

                Matrix B = B1 * B2;

                IMatrixView E = ((IContinuumMaterial3D)nurbsElement.Patch.Material).ConstitutiveMatrix;
                Matrix stiffnessMatrixGaussPoint = B.ThisTransposeTimesOtherTimesThis(E) * jacdet * gaussPoints[j].WeightFactor;

                for (int m = 0; m < nurbsElement.ControlPoints.Count * 3; m++)
                {
                    for (int n = 0; n < nurbsElement.ControlPoints.Count * 3; n++)
                    {
                        stiffnessMatrixElement[m, n] += stiffnessMatrixGaussPoint[m, n];
                    }
                }
            }
            return stiffnessMatrixElement;

        }

        private static Matrix CalculateDeformationMatrix2(IList<ControlPoint> elementControlPoints, NURBS3D nurbs, int j)
        {
            Matrix B2 = Matrix.CreateZero(9, 3 * elementControlPoints.Count);
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

        private static Matrix CalculateDeformationMatrix1(Matrix inverseJacobian)
        {
            Matrix B1 = Matrix.CreateZero(6, 9);

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

        private static Matrix CalculateInverseJacobian(Matrix jacobianMatrix, double jacdet)
        {
            Matrix inverseJacobian = Matrix.CreateZero(3, 3);

            inverseJacobian[0, 0] = jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[1, 2] * jacobianMatrix[2, 1];
            inverseJacobian[0, 1] = jacobianMatrix[0, 2] * jacobianMatrix[2, 1] - jacobianMatrix[0, 1] * jacobianMatrix[2, 2];
            inverseJacobian[0, 2] = jacobianMatrix[0, 1] * jacobianMatrix[1, 2] - jacobianMatrix[0, 2] * jacobianMatrix[1, 1];

            inverseJacobian[1, 0] = jacobianMatrix[1, 2] * jacobianMatrix[2, 0] - jacobianMatrix[1, 0] * jacobianMatrix[2, 2];
            inverseJacobian[1, 1] = jacobianMatrix[0, 0] * jacobianMatrix[2, 2] - jacobianMatrix[0, 2] * jacobianMatrix[2, 0];
            inverseJacobian[1, 2] = jacobianMatrix[0, 2] * jacobianMatrix[1, 0] - jacobianMatrix[0, 0] * jacobianMatrix[1, 2];

            inverseJacobian[2, 0] = jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[1, 1] * jacobianMatrix[2, 0];
            inverseJacobian[2, 1] = jacobianMatrix[0, 1] * jacobianMatrix[2, 0] - jacobianMatrix[0, 0] * jacobianMatrix[2, 1];
            inverseJacobian[2, 2] = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

            inverseJacobian = inverseJacobian * (1/jacdet);
            return inverseJacobian;
        }

        private static double CalculateJacobianDeterminant(Matrix jacobianMatrix)
        {
            return jacobianMatrix[0, 0] * (jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[2, 1] * jacobianMatrix[1, 2])
                                - jacobianMatrix[0, 1] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 2] - jacobianMatrix[2, 0] * jacobianMatrix[1, 2])
                                + jacobianMatrix[0, 2] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[2, 0] * jacobianMatrix[1, 1]);
        }

        private static Matrix CalculateJacobian(IList<ControlPoint> elementControlPoints, NURBS3D nurbs, int j)
        {
            Matrix jacobianMatrix = Matrix.CreateZero(3, 3);

            for (int k = 0; k < elementControlPoints.Count; k++)
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
            return gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Patch.DegreeZeta, element.Knots);
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element,Edge edge, NeumannBoundaryCondition neumann)
        {
            throw new NotSupportedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann)
        {
            throw new NotSupportedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,PressureBoundaryCondition pressure)
        {
            throw new NotSupportedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure)
        {
            throw new NotSupportedException();
        }

		public double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements)
		{
			var nurbsElement = (NURBSElement3D)element;
			var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] { element.Knots[0].Ksi, element.Knots[4].Ksi });
			var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] { element.Knots[0].Heta, element.Knots[2].Heta });
			var knotParametricCoordinatesΖeta = Vector.CreateFromArray(new double[] { element.Knots[0].Zeta, element.Knots[1].Zeta });
			NURBS3D nurbs = new NURBS3D(nurbsElement, nurbsElement.ControlPoints, knotParametricCoordinatesKsi,
				knotParametricCoordinatesHeta, knotParametricCoordinatesΖeta);
			var knotDisplacements = new double[8, 3];
			var paraviewKnotRenumbering = new int[] { 0,4,2,6,1,5,3,7 };
			for (int j = 0; j < element.Knots.Count; j++)
			{
				for (int i = 0; i < element.ControlPoints.Count; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
					knotDisplacements[paraviewKnotRenumbering[j], 1] += nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
					knotDisplacements[paraviewKnotRenumbering[j], 2] += nurbs.NurbsValues[i, j] * localDisplacements[i, 2];
				}
			}

			return knotDisplacements;
		}
	}
}
