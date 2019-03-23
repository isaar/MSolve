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
    public class NURBS3D
    {
        public IMatrix2D NurbsValues { get; private set; }
        public IMatrix2D NurbsDerivativeValuesKsi { get; private set; }
        public IMatrix2D NurbsDerivativeValuesHeta { get; private set; }
        public IMatrix2D NurbsDerivativeValuesZeta { get; private set; }

        public NURBS3D(Element element, IList<ControlPoint> controlPoints)
        {
            GaussQuadrature gauss = new GaussQuadrature();
            IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Patch.DegreeZeta, element.Knots);

            IVector parametricGaussPointKsi = new Vector(element.Patch.DegreeKsi + 1);
            for (int i = 0; i < element.Patch.DegreeKsi + 1; i++)
            {
                parametricGaussPointKsi[i] = gaussPoints[i * (element.Patch.DegreeZeta + 1) * (element.Patch.DegreeHeta + 1)].Ksi;
            }

            IVector parametricGaussPointHeta = new Vector(element.Patch.DegreeHeta + 1);
            for (int i = 0; i < element.Patch.DegreeHeta + 1; i++)
            {
                parametricGaussPointHeta[i] = gaussPoints[i * (element.Patch.DegreeZeta + 1)].Heta;
            }

            IVector parametricGaussPointZeta = new Vector(element.Patch.DegreeZeta + 1);
            for (int i = 0; i < element.Patch.DegreeZeta + 1; i++)
            {
                parametricGaussPointZeta[i] = gaussPoints[i].Zeta;
            }

            BSPLines1D bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi, parametricGaussPointKsi);
            BSPLines1D bsplinesHeta = new BSPLines1D(element.Patch.DegreeHeta, element.Patch.KnotValueVectorHeta, parametricGaussPointHeta);
            BSPLines1D bsplinesZeta = new BSPLines1D(element.Patch.DegreeZeta, element.Patch.KnotValueVectorZeta, parametricGaussPointZeta);
            bsplinesKsi.calculateBSPLinesAndDerivatives();
            bsplinesHeta.calculateBSPLinesAndDerivatives();
            bsplinesZeta.calculateBSPLinesAndDerivatives();

            int supportKsi = element.Patch.DegreeKsi + 1;
            int supportHeta = element.Patch.DegreeHeta + 1;
            int supportZeta = element.Patch.DegreeZeta + 1;
            int numberOfElementControlPoints = supportKsi * supportHeta * supportZeta;

            NurbsValues = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
            NurbsDerivativeValuesKsi = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
            NurbsDerivativeValuesHeta = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
            NurbsDerivativeValuesZeta = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);

            for (int i = 0; i < supportKsi; i++)
            {
                for (int j = 0; j < supportHeta; j++)
                {
                    for (int k = 0; k < supportZeta; k++)
                    {
                        double sumKsiHetaZeta = 0;
                        double sumdKsiHetaZeta = 0;
                        double sumKsidHetaZeta = 0;
                        double sumKsiHetadZeta = 0;

                        for (int m = 0; m < numberOfElementControlPoints; m++)
                        {
                            int indexKsi = element.ControlPoints[m].ID / (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta);
                            int indexHeta = element.ControlPoints[m].ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) / element.Patch.NumberOfControlPointsZeta;
                            int indexZeta = element.ControlPoints[m].ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) % element.Patch.NumberOfControlPointsZeta;

                            sumKsiHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;

                            sumdKsiHetaZeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;

                            sumKsidHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;

                            sumKsiHetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;
                        }
                        for (int m = 0; m < numberOfElementControlPoints; m++)
                        {
                            int indexKsi = element.ControlPoints[m].ID / (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta);
                            int indexHeta = element.ControlPoints[m].ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) / element.Patch.NumberOfControlPointsZeta;
                            int indexZeta = element.ControlPoints[m].ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) % element.Patch.NumberOfControlPointsZeta;

                            NurbsValues[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / sumKsiHetaZeta;

                            NurbsDerivativeValuesKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHetaZeta -
                                bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHetaZeta) *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

                            NurbsDerivativeValuesHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                (bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHetaZeta -
                                bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHetaZeta) *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

                            NurbsDerivativeValuesZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                (bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] * sumKsiHetaZeta -
                                bsplinesZeta.BSPLineValues[indexZeta, k] * sumKsiHetadZeta) *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);
                        }
                    }
                }
            }
        }

	    public NURBS3D(Element element, IList<ControlPoint> controlPoints, Vector parametricGaussPointKsi,
		    Vector parametricGaussPointHeta, Vector parametricGaussPointZeta)
	    {
		    var parametricPointsCount = parametricGaussPointKsi.Length * parametricGaussPointHeta.Length *
		                                parametricGaussPointZeta.Length;

			BSPLines1D bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi,
			    parametricGaussPointKsi);
		    BSPLines1D bsplinesHeta = new BSPLines1D(element.Patch.DegreeHeta, element.Patch.KnotValueVectorHeta,
			    parametricGaussPointHeta);
		    BSPLines1D bsplinesZeta = new BSPLines1D(element.Patch.DegreeZeta, element.Patch.KnotValueVectorZeta,
			    parametricGaussPointZeta);
		    bsplinesKsi.calculateBSPLinesAndDerivatives();
		    bsplinesHeta.calculateBSPLinesAndDerivatives();
		    bsplinesZeta.calculateBSPLinesAndDerivatives();

		    int supportKsi = parametricGaussPointKsi.Length;
		    int supportHeta = parametricGaussPointHeta.Length;
		    int supportZeta = parametricGaussPointZeta.Length;
		    int numberOfElementControlPoints = (element.Patch.DegreeKsi + 1) * (element.Patch.DegreeHeta + 1) *
		                                       (element.Patch.DegreeZeta + 1);

		    NurbsValues = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);
		    NurbsDerivativeValuesKsi = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);
		    NurbsDerivativeValuesHeta = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);
		    NurbsDerivativeValuesZeta = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);

		    for (int i = 0; i < supportKsi; i++)
		    {
			    for (int j = 0; j < supportHeta; j++)
			    {
				    for (int k = 0; k < supportZeta; k++)
				    {
					    double sumKsiHetaZeta = 0;
					    double sumdKsiHetaZeta = 0;
					    double sumKsidHetaZeta = 0;
					    double sumKsiHetadZeta = 0;

					    for (int m = 0; m < numberOfElementControlPoints; m++)
					    {
						    int indexKsi = element.ControlPoints[m].ID /
						                   (element.Patch.NumberOfControlPointsHeta *
						                    element.Patch.NumberOfControlPointsZeta);
						    int indexHeta = element.ControlPoints[m].ID %
						                    (element.Patch.NumberOfControlPointsHeta *
						                     element.Patch.NumberOfControlPointsZeta) /
						                    element.Patch.NumberOfControlPointsZeta;
						    int indexZeta = element.ControlPoints[m].ID %
						                    (element.Patch.NumberOfControlPointsHeta *
						                     element.Patch.NumberOfControlPointsZeta) %
						                    element.Patch.NumberOfControlPointsZeta;

						    sumKsiHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
						                      bsplinesHeta.BSPLineValues[indexHeta, j] *
						                      bsplinesZeta.BSPLineValues[indexZeta, k] *
						                      controlPoints[m].WeightFactor;

						    sumdKsiHetaZeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
						                       bsplinesHeta.BSPLineValues[indexHeta, j] *
						                       bsplinesZeta.BSPLineValues[indexZeta, k] *
						                       controlPoints[m].WeightFactor;

						    sumKsidHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
						                       bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
						                       bsplinesZeta.BSPLineValues[indexZeta, k] *
						                       controlPoints[m].WeightFactor;

						    sumKsiHetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
						                       bsplinesHeta.BSPLineValues[indexHeta, j] *
						                       bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
						                       controlPoints[m].WeightFactor;
					    }

					    for (int m = 0; m < numberOfElementControlPoints; m++)
					    {
						    int indexKsi = element.ControlPoints[m].ID /
						                   (element.Patch.NumberOfControlPointsHeta *
						                    element.Patch.NumberOfControlPointsZeta);
						    int indexHeta = element.ControlPoints[m].ID %
						                    (element.Patch.NumberOfControlPointsHeta *
						                     element.Patch.NumberOfControlPointsZeta) /
						                    element.Patch.NumberOfControlPointsZeta;
						    int indexZeta = element.ControlPoints[m].ID %
						                    (element.Patch.NumberOfControlPointsHeta *
						                     element.Patch.NumberOfControlPointsZeta) %
						                    element.Patch.NumberOfControlPointsZeta;

						    NurbsValues[m, i * supportHeta * supportZeta + j * supportZeta + k] =
							    bsplinesKsi.BSPLineValues[indexKsi, i] *
							    bsplinesHeta.BSPLineValues[indexHeta, j] *
							    bsplinesZeta.BSPLineValues[indexZeta, k] *
							    controlPoints[m].WeightFactor / sumKsiHetaZeta;

						    NurbsDerivativeValuesKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
							    (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHetaZeta -
							     bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHetaZeta) *
							    bsplinesHeta.BSPLineValues[indexHeta, j] *
							    bsplinesZeta.BSPLineValues[indexZeta, k] *
							    controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

						    NurbsDerivativeValuesHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
							    bsplinesKsi.BSPLineValues[indexKsi, i] *
							    (bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHetaZeta -
							     bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHetaZeta) *
							    bsplinesZeta.BSPLineValues[indexZeta, k] *
							    controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

						    NurbsDerivativeValuesZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
							    bsplinesKsi.BSPLineValues[indexKsi, i] *
							    bsplinesHeta.BSPLineValues[indexHeta, j] *
							    (bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] * sumKsiHetaZeta -
							     bsplinesZeta.BSPLineValues[indexZeta, k] * sumKsiHetadZeta) *
							    controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);
					    }
				    }
			    }
		    }
	    }
    }
}
