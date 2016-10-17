using System;
using System.Collections.Generic;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.PreProcessor.Elements.SupportiveClasses;
using ISAAR.MSolve.Matrices;
using ISSAR.MSolve.IGAPreProcessor.Integration;
using ISSAR.MSolve.IGAPreProcessor.Interfaces;

namespace ISSAR.MSolve.IGAPreProcessor.Elements
{
    public class NURBSElement2D:IIsogeometricStructuralElement
    {
        private int id;
        private IList<Knot> knots;
        public IVector<int> connectivity { get; }
        private IGAModel model;
        public IList<GaussLegendrePoint3D> gaussPoints { get; }
        
        public NURBSElement2D(int id,IGAModel model, IList<Knot> knots, IVector<int> connectivity)
        {
            this.id = id;
            this.knots = knots;
            this.connectivity = connectivity;
            this.model = model;

            this.CreateElementGaussPoints();
        }
        #region IISogeometricStructuralElement
        public int ID
        {
            get
            {
                return this.id;
            }
        }

        public ElementDimensions ElementDimensions
        {
            get
            {
                return ElementDimensions.TwoD;
            }
        }

        public IMatrix2D<double> StiffnessMatrix(IGAElement element)
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(IGAElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(IGAElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(IGAElement element, double[] localDisplacements)
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

            IVector<double> coordinatesKsi = new Vector<double>(gaussPointsPerAxisKsi.Length);
            IVector<double> weightsKsi = new Vector<double>(gaussPointsPerAxisKsi.Length);
            for (int indexKsi = 0; indexKsi < gaussPointsPerAxisKsi.Length; indexKsi++)
            {
                coordinatesKsi[indexKsi] = 0.5 * (knots[0].Ksi + knots[2].Ksi+ (knots[2].Ksi - knots[0].Ksi) * gaussPointsPerAxisKsi[indexKsi].Coordinate);
                weightsKsi[indexKsi] = 0.5* ((knots[2].Ksi - knots[0].Ksi) * gaussPointsPerAxisKsi[indexKsi].WeightFactor);
            }

            IVector<double> coordinatesHeta = new Vector<double>(gaussPointsPerAxisHeta.Length);
            IVector<double> weightsHeta = new Vector<double>(gaussPointsPerAxisHeta.Length);
            for (int indexHeta = 0; indexHeta < gaussPointsPerAxisHeta.Length; indexHeta++)
            {
                coordinatesHeta[indexHeta] = 0.5 * (knots[0].Heta + knots[2].Heta + (knots[2].Heta - knots[0].Heta) * gaussPointsPerAxisHeta[indexHeta].Coordinate);
                weightsHeta[indexHeta] = 0.5 * ((knots[2].Heta - knots[0].Heta) * gaussPointsPerAxisHeta[indexHeta].WeightFactor);
            }

            IList<ControlPoint> controlPointsElement = new List<ControlPoint>();
            for (int k = 0; k < (model.DegreeKsi+1)*(model.DegreeHeta+1); k++)
            {
                controlPointsElement.Add(model.ControlPoints[this.connectivity[k]]);
            }

            ShapeNURBS2D nurbs = new ShapeNURBS2D(this, model, coordinatesKsi, coordinatesHeta, controlPointsElement);

            for (int i = 0; i < gaussPointsPerAxisKsi.Length; i++)
            {
                for (int j = 0; j < gaussPointsPerAxisHeta.Length; j++)
                {
                    IMatrix2D<double> jacobianMatrix = new Matrix2D<double>(2, 2);

                    for (int k = 0; k < controlPointsElement.Count; k++)
                    {
                        jacobianMatrix[0, 0] += nurbs.nurbsDerivativeValuesKsi[k, j] * controlPointsElement[k].X;
                        jacobianMatrix[0, 1] += nurbs.nurbsDerivativeValuesKsi[k, j] * controlPointsElement[k].Y;
                        jacobianMatrix[1, 0] += nurbs.nurbsDerivativeValuesHeta[k, j] * controlPointsElement[k].X;
                        jacobianMatrix[1, 1] += nurbs.nurbsDerivativeValuesHeta[k, j] * controlPointsElement[k].Y;
                    }

                    double jacdet = jacobianMatrix[0, 0] * jacobianMatrix[1, 1]
                    - jacobianMatrix[1, 0] * jacobianMatrix[0, 1];


                    Matrix2D<double> B1 = new Matrix2D<double>(3, 4);

                    B1[0, 0] += jacobianMatrix[1, 1] / jacdet;
                    B1[0, 1] += -jacobianMatrix[0, 1] / jacdet;
                    B1[1, 2] += -jacobianMatrix[1, 0] / jacdet;
                    B1[1, 3] += jacobianMatrix[0, 0] / jacdet;
                    B1[2, 0] += -jacobianMatrix[1, 0] / jacdet;
                    B1[2, 1] += jacobianMatrix[0, 0] / jacdet;
                    B1[2, 2] += jacobianMatrix[1, 1] / jacdet;
                    B1[2, 3] += -jacobianMatrix[0, 1] / jacdet;

                    Matrix2D<double> B2 = new Matrix2D<double>(4, 2 * controlPointsElement.Count);

                    for (int column = 0; column < 2*controlPointsElement.Count; column+=2)
                    {
                        B2[0, column] += nurbs.nurbsDerivativeValuesKsi[column / 2, j];
                        B2[1, column] += nurbs.nurbsDerivativeValuesHeta[column / 2, j];
                        B2[2, column + 1] += nurbs.nurbsDerivativeValuesKsi[column / 2, j];
                        B2[3, column + 1] += nurbs.nurbsDerivativeValuesHeta[column / 2, j];
                    }

                    IMatrix2D<double> B = B1 * B2;


                    double weightFactor = gaussPointsPerAxisKsi[i].WeightFactor * gaussPointsPerAxisHeta[i].WeightFactor * jacdet;
                    this.gaussPoints.Add(new GaussLegendrePoint3D(gaussPointsPerAxisKsi[i].Coordinate, gaussPointsPerAxisHeta[j].Coordinate, 0.0,B, weightFactor));
                }
            }
            



        }

    }
}
