using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Elements
{
    public class NURBSElement2DCollocationBoundary:Element, IStructuralIsogeometricElement, ICollocationElement
    {
        protected readonly static IDofType[] controlPointDOFTypes = new StructuralDof[] { StructuralDof.TranslationX, StructuralDof.TranslationY };
        protected IDofType[][] dofTypes;
        private CollocationPoint2D _collocationPoint;

        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

        public IElementDofEnumerator DofEnumerator
        {
            get { return dofEnumerator; }
            set { this.dofEnumerator = value; }
        }

        public bool MaterialModified { get; }

        INode ICollocationElement.CollocationPoint { get => _collocationPoint; set => _collocationPoint = (CollocationPoint2D)value; }

        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
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
            var elementCollocation = (NURBSElement2DCollocation)element;

            var nurbs = new NURBS2D(elementCollocation.Patch.DegreeKsi, elementCollocation.Patch.DegreeHeta,
                elementCollocation.Patch.KnotValueVectorKsi, elementCollocation.Patch.KnotValueVectorHeta,
                elementCollocation.CollocationPoint, elementCollocation.ControlPoints);

            var jacobianMatrix = CalculateJacobianMatrix(elementCollocation, nurbs);

            var (xGaussPoint, yGaussPoint) = CalculateNormalVectors(elementCollocation, nurbs);

            var inverseJacobian = jacobianMatrix.Invert();
            var dR = CalculateNaturalDerivatives(nurbs, inverseJacobian);

            return CalculateCollocationPointStiffness(elementCollocation, xGaussPoint, dR, yGaussPoint);
            
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

        private Matrix CalculateCollocationPointStiffness(NURBSElement2DCollocation elementCollocation,
            double xGaussPoint, double[,] dR, double yGaussPoint)
        {
            var collocationPointStiffness = Matrix.CreateZero(2, elementCollocation.ControlPoints.Count * 2);
            var E = elementCollocation.Patch.Material.YoungModulus;
            var nu = elementCollocation.Patch.Material.PoissonRatio;
            var lambda = E * nu / (1 + nu) / (1 - 2 * nu);
            var m = E / (2 * (1 + nu));

            for (int i = 0; i < elementCollocation.ControlPoints.Count * 2; i += 2)
            {
                var index = i / 2;
                collocationPointStiffness[0, i] =
                    (lambda + 2 * m) * xGaussPoint * dR[0, index] + m * yGaussPoint * dR[1, index];
                collocationPointStiffness[1, i] = lambda * xGaussPoint * dR[1, index] + m * yGaussPoint * dR[0, index];

                collocationPointStiffness[0, i + 1] = lambda * yGaussPoint * dR[0, index] + m * xGaussPoint * dR[1, index];
                collocationPointStiffness[1, i + 1] =
                    (lambda + 2 * m) * yGaussPoint * dR[1, index] + m * xGaussPoint * dR[0, index];
            }

            return collocationPointStiffness;
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
        
        public IMatrix MassMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public IMatrix DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public IList<IList<IDofType>> GetElementDOFTypes(IElement element)
        {
            throw new NotImplementedException();
        }

        public INode CollocationPoint { get; set; }
    }
}
