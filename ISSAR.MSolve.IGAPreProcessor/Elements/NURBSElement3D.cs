using ISSAR.MSolve.IGAPreProcessor.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.PreProcessor.Elements.SupportiveClasses;
using ISAAR.MSolve.Matrices;
using ISSAR.MSolve.IGAPreProcessor.Integration;

namespace ISSAR.MSolve.IGAPreProcessor.Elements
{
    class NURBSElement3D : IIsogeometricStructuralElement
    {
        private int id;
        private IList<Knot> knots;
        public IVector<int> connectivity { get; }
        private IGAModel model;
        public IList<GaussLegendrePoint3D> gaussPoints { get; }

        public NURBSElement3D(int id, IGAModel model, IList<Knot> knots, IVector<int> connectivity)
        {
            this.id = id;
            this.knots = knots;
            this.connectivity = connectivity;
            this.model = model;

            this.CreateElementGaussPoints();
        }

        #region IISogeometricStructuralElement
        public ElementDimensions ElementDimensions
        {
            get
            {
               return ElementDimensions.ThreeD;
            }
        }

        public int ID
        {
            get
            {
                return this.id;
            }
        }

        public double[] CalculateForces(IGAElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(IGAElement element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(IGAElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public IMatrix2D<double> StiffnessMatrix(IGAElement element)
        {
            throw new NotImplementedException();
        }
        #endregion

        private void CreateElementGaussPoints()
        {
            GaussLegendrePoint1D[] gaussPointsPerAxisKsi =
                GaussQuadrature.GetGaussLegendrePoints(model.DegreeKsi);

            GaussLegendrePoint1D[] gaussPointsPerAxisHeta =
                GaussQuadrature.GetGaussLegendrePoints(model.DegreeHeta);

            GaussLegendrePoint1D[] gaussPointsPerAxisZeta =
                GaussQuadrature.GetGaussLegendrePoints(model.DegreeZeta);

            IVector<double> coordinatesKsi = new Vector<double>(gaussPointsPerAxisKsi.Length);
            IVector<double> weightsKsi = new Vector<double>(gaussPointsPerAxisKsi.Length);
            for (int indexKsi = 0; indexKsi < gaussPointsPerAxisKsi.Length; indexKsi++)
            {
                coordinatesKsi[indexKsi] = 0.5 * (knots[0].Ksi + knots[4].Ksi + (knots[4].Ksi - knots[0].Ksi) * gaussPointsPerAxisKsi[indexKsi].Coordinate);
                weightsKsi[indexKsi] = 0.5 * ((knots[4].Ksi - knots[0].Ksi) * gaussPointsPerAxisKsi[indexKsi].WeightFactor);
            }

            IVector<double> coordinatesHeta = new Vector<double>(gaussPointsPerAxisHeta.Length);
            IVector<double> weightsHeta = new Vector<double>(gaussPointsPerAxisHeta.Length);
            for (int indexHeta = 0; indexHeta < gaussPointsPerAxisHeta.Length; indexHeta++)
            {
                coordinatesHeta[indexHeta] = 0.5 * (knots[0].Heta + knots[2].Heta + (knots[2].Heta - knots[0].Heta) * gaussPointsPerAxisHeta[indexHeta].Coordinate);
                weightsHeta[indexHeta] = 0.5 * ((knots[2].Heta - knots[0].Heta) * gaussPointsPerAxisHeta[indexHeta].WeightFactor);
            }

            IVector<double> coordinatesZeta = new Vector<double>(gaussPointsPerAxisZeta.Length);
            IVector<double> weightsZeta = new Vector<double>(gaussPointsPerAxisZeta.Length);
            for (int indexZeta = 0; indexZeta < gaussPointsPerAxisZeta.Length; indexZeta++)
            {
                coordinatesZeta[indexZeta] = 0.5 * (knots[0].Zeta + knots[1].Zeta + (knots[1].Zeta - knots[0].Zeta) * gaussPointsPerAxisZeta[indexZeta].Coordinate);
                weightsZeta[indexZeta] = 0.5 * ((knots[1].Zeta - knots[0].Zeta) * gaussPointsPerAxisZeta[indexZeta].WeightFactor);
            }

            IList<ControlPoint> controlPointsElement = new List<ControlPoint>();
            for (int k = 0; k < (model.DegreeKsi + 1) * (model.DegreeHeta + 1) * (model.DegreeZeta + 1); k++)
            {
                controlPointsElement.Add(model.ControlPoints[this.connectivity[k]]);
            }

            ShapeNURBS3D nurbs = new ShapeNURBS3D(this, model, coordinatesKsi,coordinatesHeta,coordinatesZeta,controlPointsElement);

            for (int j = 0; j < gaussPointsPerAxisKsi.Length*gaussPointsPerAxisHeta.Length*gaussPointsPerAxisZeta.Length; j++)
            {
                Matrix2D<double> jacobianMatrix = new Matrix2D<double>(3, 3);

                for (int k = 0; k < controlPointsElement.Count; k++)
                {
                    jacobianMatrix[0, 0] += nurbs.nurbsDerivativeValuesKsi[k, j] * controlPointsElement[k].X;
                    jacobianMatrix[0, 1] += nurbs.nurbsDerivativeValuesKsi[k, j] * controlPointsElement[k].Y;
                    jacobianMatrix[0, 2] += nurbs.nurbsDerivativeValuesKsi[k, j] * controlPointsElement[k].Z;
                    jacobianMatrix[1, 0] += nurbs.nurbsDerivativeValuesHeta[k, j] * controlPointsElement[k].X;
                    jacobianMatrix[1, 1] += nurbs.nurbsDerivativeValuesHeta[k, j] * controlPointsElement[k].Y;
                    jacobianMatrix[1, 2] += nurbs.nurbsDerivativeValuesHeta[k, j] * controlPointsElement[k].Z;
                    jacobianMatrix[2, 0] += nurbs.nurbsDerivativeValuesZeta[k, j] * controlPointsElement[k].X;
                    jacobianMatrix[2, 1] += nurbs.nurbsDerivativeValuesZeta[k, j] * controlPointsElement[k].Y;
                    jacobianMatrix[2, 2] += nurbs.nurbsDerivativeValuesZeta[k, j] * controlPointsElement[k].Z;

                }

                double jacdet = jacobianMatrix[0, 0] * (jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[2, 1] * jacobianMatrix[1, 2])
                    - jacobianMatrix[0, 1] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 2] - jacobianMatrix[2, 0] * jacobianMatrix[1, 2])
                    + jacobianMatrix[0, 2] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[2, 0] * jacobianMatrix[1, 1]);

                Matrix2D<double> inverseJacobian = new Matrix2D<double>(3, 3);

                inverseJacobian[0, 0] = jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[1, 2] * jacobianMatrix[2, 1];
                inverseJacobian[0, 1] = jacobianMatrix[0, 2] * jacobianMatrix[2, 1] - jacobianMatrix[0, 1] * jacobianMatrix[2, 2];
                inverseJacobian[0, 2] = jacobianMatrix[0, 1] * jacobianMatrix[1, 2] - jacobianMatrix[0, 2] * jacobianMatrix[1, 1];
                inverseJacobian[1, 0] = jacobianMatrix[1, 2] * jacobianMatrix[2, 0] - jacobianMatrix[1, 0] * jacobianMatrix[2, 2];
                inverseJacobian[1, 1] = jacobianMatrix[0, 0] * jacobianMatrix[2, 2] - jacobianMatrix[0, 2] * jacobianMatrix[2, 0];
                inverseJacobian[1, 2] = jacobianMatrix[0, 2] * jacobianMatrix[1, 0] - jacobianMatrix[0, 0] * jacobianMatrix[1, 2];
                inverseJacobian[2, 0] = jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[1, 1] * jacobianMatrix[2, 0];
                inverseJacobian[2, 1] = jacobianMatrix[0, 1] * jacobianMatrix[2, 0] - jacobianMatrix[0, 0] * jacobianMatrix[2, 1];
                inverseJacobian[2, 2] = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];
                inverseJacobian =inverseJacobian * (1 / jacdet);

                Matrix2D<double> B1 = new Matrix2D<double>(6, 9);

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

                Matrix2D<double> B2 = new Matrix2D<double>(9, 3 * controlPointsElement.Count);
                for (int column = 0; column < 3 * controlPointsElement.Count; column += 3)
                {
                    B2[0, column] += nurbs.nurbsDerivativeValuesKsi[column / 3, j];
                    B2[1, column] += nurbs.nurbsDerivativeValuesHeta[column / 3, j];
                    B2[2, column] += nurbs.nurbsDerivativeValuesZeta[column / 3, j];

                    B2[3, column + 1] += nurbs.nurbsDerivativeValuesKsi[column / 3, j];
                    B2[4, column + 1] += nurbs.nurbsDerivativeValuesHeta[column / 3, j];
                    B2[5, column + 1] += nurbs.nurbsDerivativeValuesZeta[column / 3, j];

                    B2[6, column + 2] += nurbs.nurbsDerivativeValuesKsi[column / 3, j];
                    B2[7, column + 2] += nurbs.nurbsDerivativeValuesHeta[column / 3, j];
                    B2[8, column + 2] += nurbs.nurbsDerivativeValuesZeta[column / 3, j];
                }

                Matrix2D<double> B = B1 * B2;

                int indexKsi = j / (gaussPointsPerAxisHeta.Length * gaussPointsPerAxisZeta.Length);
                int indexHeta = j % (gaussPointsPerAxisHeta.Length * gaussPointsPerAxisZeta.Length) / gaussPointsPerAxisZeta.Length;
                int indexZeta = j % (gaussPointsPerAxisHeta.Length * gaussPointsPerAxisZeta.Length) % gaussPointsPerAxisZeta.Length;

                double weightFactor = gaussPointsPerAxisKsi[indexKsi].WeightFactor *
                    gaussPointsPerAxisHeta[indexHeta].WeightFactor *
                    gaussPointsPerAxisZeta[indexZeta].WeightFactor * jacdet;
                this.gaussPoints.Add(new GaussLegendrePoint3D(gaussPointsPerAxisKsi[indexKsi].Coordinate, gaussPointsPerAxisHeta[indexHeta].Coordinate, gaussPointsPerAxisZeta[indexZeta].Coordinate, B, weightFactor));

            }











        }
    }
}
