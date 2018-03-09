using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Problems.Structural.Elements;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.IGA.Problems.SupportiveClasses
{
    public class NURBS2D
    {
        public IMatrix2D NurbsValues { get; private set; }
        public IMatrix2D NurbsDerivativeValuesKsi { get; private set; }
        public IMatrix2D NurbsDerivativeValuesHeta { get; private set; }

        public NURBS2D(Element element, IList<ControlPoint> controlPoints)
        {
            GaussQuadrature gauss = new GaussQuadrature();
            IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Knots);

            IVector parametricGaussPointKsi = new Vector(element.Patch.DegreeKsi + 1);
            for (int i = 0; i < element.Patch.DegreeKsi+1; i++)
            {
                parametricGaussPointKsi[i] = gaussPoints[i * (element.Patch.DegreeHeta + 1)].Ksi;
            }

            IVector parametricGaussPointHeta = new Vector(element.Patch.DegreeHeta + 1);
            for (int i = 0; i < element.Patch.DegreeHeta + 1; i++)
            {
                parametricGaussPointHeta[i] = gaussPoints[i ].Heta;
            }

            BSPLines1D bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi, parametricGaussPointKsi);
            BSPLines1D bsplinesHeta = new BSPLines1D(element.Patch.DegreeHeta, element.Patch.KnotValueVectorHeta, parametricGaussPointHeta);
            bsplinesKsi.calculateBSPLinesAndDerivatives();
            bsplinesHeta.calculateBSPLinesAndDerivatives();
            

            int supportKsi = element.Patch.DegreeKsi + 1;
            int supportHeta = element.Patch.DegreeHeta + 1;
            int numberOfElementControlPoints = supportKsi * supportHeta;

            NurbsValues = new Matrix2D(numberOfElementControlPoints,gaussPoints.Count);
            NurbsDerivativeValuesKsi = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
            NurbsDerivativeValuesHeta = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);

            for (int i = 0; i < supportKsi; i++)
            {
                for (int j = 0; j < supportHeta; j++)
                {
                    double sumKsiHeta = 0;
                    double sumdKsiHeta = 0;
                    double sumKsidHeta = 0;

                    for (int k = 0; k < numberOfElementControlPoints; k++)
                    {
                        // Why type casting is needed.?

                        int indexKsi = element.ControlPoints[k].ID / element.Patch.NumberOfControlPointsHeta;
                        int indexHeta = element.ControlPoints[k].ID % element.Patch.NumberOfControlPointsHeta;
                        sumKsiHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                            bsplinesHeta.BSPLineValues[indexHeta, j] *
                            controlPoints[k].WeightFactor;
                        sumdKsiHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i]*
                            bsplinesHeta.BSPLineValues[indexHeta, j] *
                            controlPoints[k].WeightFactor;
                        sumKsidHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                            bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                            controlPoints[k].WeightFactor;
                    }
                    for (int k = 0; k < numberOfElementControlPoints; k++)
                    {
                        int indexKsi = element.ControlPoints[k].ID / element.Patch.NumberOfControlPointsHeta;
                        int indexHeta = element.ControlPoints[k].ID % element.Patch.NumberOfControlPointsHeta;

                        NurbsValues[k, i * supportHeta + j] = 
                            bsplinesKsi.BSPLineValues[indexKsi, i] *
                            bsplinesHeta.BSPLineValues[indexHeta, j] *
                            controlPoints[k].WeightFactor / sumKsiHeta;

                        NurbsDerivativeValuesKsi[k, i * supportHeta + j] =
                            bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
                            (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHeta -
                            bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHeta) / Math.Pow(sumKsiHeta, 2);

                        NurbsDerivativeValuesHeta[k, i * supportHeta + j] =
                            bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
                            (bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHeta -
                            bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHeta) / Math.Pow(sumKsiHeta, 2);
                    }
                }
            }
        }
    }
}
