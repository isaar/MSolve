using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.IGA.Problems.SupportiveClasses
{
    public class NURBS1D
    {
        public IMatrix2D NurbsValues { get; private set; }
        public IMatrix2D NurbsDerivativeValuesKsi { get; private set; }

        public NURBS1D(Element element, IList<ControlPoint> controlPoints)
        {
            GaussQuadrature gauss = new GaussQuadrature();
            IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Knots);

            IVector parametricGaussPointKsi = new Vector(element.Patch.DegreeKsi + 1);
            for (int i = 0; i < element.Patch.DegreeKsi + 1; i++)
            {
                parametricGaussPointKsi[i] = gaussPoints[i].Ksi;
            }
            BSPLines1D bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi, parametricGaussPointKsi);
            bsplinesKsi.calculateBSPLinesAndDerivatives();

            int supportKsi = element.Patch.DegreeKsi + 1;
            int numberOfElementControlPoints = supportKsi;

            NurbsValues = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
            NurbsDerivativeValuesKsi = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
            for (int i = 0; i < supportKsi; i++)
            {
                double sumKsi = 0;
                double sumdKsi = 0;

                for (int j = 0; j < numberOfElementControlPoints; j++)
                {
                    int indexKsi = element.ControlPoints[j].ID;
                    sumKsi += bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[j].WeightFactor;
                    sumdKsi += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * controlPoints[j].WeightFactor;
                }
                for (int j = 0; j < numberOfElementControlPoints; j++)
                {
                    int indexKsi = element.ControlPoints[j].ID;
                    NurbsValues[j, i] = bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[j].WeightFactor / sumKsi;
                    NurbsDerivativeValuesKsi[j, i] = controlPoints[j].WeightFactor*(bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsi -
                        bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsi) / Math.Pow(sumKsi, 2);
                }
            }
        }


        public NURBS1D(Element element, IList<ControlPoint> controlPoints, Edge edge)
        {
            GaussQuadrature gauss = new GaussQuadrature();
            IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(edge.Degree, element.Knots);

            IVector parametricGaussPointKsi = new Vector(edge.Degree + 1);
            for (int i = 0; i < edge.Degree + 1; i++)
            {
                parametricGaussPointKsi[i] = gaussPoints[i].Ksi;
            }
            BSPLines1D bsplinesKsi = new BSPLines1D(edge.Degree, edge.KnotValueVector, parametricGaussPointKsi);
            bsplinesKsi.calculateBSPLinesAndDerivatives();

            int supportKsi = edge.Degree + 1;
            int numberOfElementControlPoints = supportKsi;

            NurbsValues = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
            NurbsDerivativeValuesKsi = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);

            for (int i = 0; i < supportKsi; i++)
            {
                double sumKsi = 0;
                double sumdKsi = 0;

                for (int j = 0; j < numberOfElementControlPoints; j++)
                {
                    int index = controlPoints[j].ID;
                    sumKsi += bsplinesKsi.BSPLineValues[index, i] * controlPoints[j].WeightFactor;
                    sumdKsi += bsplinesKsi.BSPLineDerivativeValues[index, i] * controlPoints[j].WeightFactor;
                }
                for (int j = 0; j < numberOfElementControlPoints; j++)
                {
                    int indexKsi = controlPoints[j].ID;
                    NurbsValues[j, i] = bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[j].WeightFactor / sumKsi;
                    NurbsDerivativeValuesKsi[j, i] = controlPoints[j].WeightFactor * (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsi -
                        bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsi) / Math.Pow(sumKsi, 2);
                }
            }
        }


    }
}
