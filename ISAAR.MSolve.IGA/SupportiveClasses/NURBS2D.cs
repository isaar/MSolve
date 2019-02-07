using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.IGA.Elements;

namespace ISAAR.MSolve.IGA.Problems.SupportiveClasses
{
    public class NURBS2D
    {
		public IMatrix2D NurbsValues { get; private set; }
        public IMatrix2D NurbsDerivativeValuesKsi { get; private set; }
        public IMatrix2D NurbsDerivativeValuesHeta { get; private set; }
	    public IMatrix2D NurbsSecondDerivativeValueKsi { get; private set; }
		public IMatrix2D NurbsSecondDerivativeValueHeta { get; private set; }
		public IMatrix2D NurbsSecondDerivativeValueKsiHeta { get; private set; }

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
	        NurbsSecondDerivativeValueKsi = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
			NurbsSecondDerivativeValueHeta = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
			NurbsSecondDerivativeValueKsiHeta = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);

			for (int i = 0; i < supportKsi; i++)
            {
                for (int j = 0; j < supportHeta; j++)
                {
                    double sumKsiHeta = 0;
                    double sumdKsiHeta = 0;
                    double sumKsidHeta = 0;
	                double sumdKsidKsi = 0;
	                double sumdHetadHeta = 0;
	                double sumdKsidHeta = 0;

					for (int k = 0; k < numberOfElementControlPoints; k++)
                    {
                        // Why type casting is needed.?

                        int indexKsi = element.ControlPoints[k].ID / element.Patch.NumberOfControlPointsHeta;
                        int indexHeta = element.ControlPoints[k].ID % element.Patch.NumberOfControlPointsHeta;
	                    sumKsiHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
	                                  bsplinesHeta.BSPLineValues[indexHeta, j] *
	                                  controlPoints[k].WeightFactor;
	                    sumdKsiHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
	                                   bsplinesHeta.BSPLineValues[indexHeta, j] *
	                                   controlPoints[k].WeightFactor;
	                    sumKsidHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
	                                   bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
	                                   controlPoints[k].WeightFactor;
	                    sumdKsidKsi += bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] *
	                                   bsplinesHeta.BSPLineValues[indexHeta, j] *
	                                   controlPoints[k].WeightFactor;
	                    sumdHetadHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
	                                     bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] *
	                                     controlPoints[k].WeightFactor;
	                    sumdKsidHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
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

	                    NurbsSecondDerivativeValueKsi[k, i * supportHeta + j] =
		                    bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
		                    (bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] / sumKsiHeta -
		                     2 * bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
		                     bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
		                     2 * bsplinesKsi.BSPLineValues[indexKsi, i] * Math.Pow(sumdKsiHeta, 2) / Math.Pow(sumKsiHeta, 3));

	                    NurbsSecondDerivativeValueHeta[k, i * supportHeta + j] =
		                    bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
		                    (bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] / sumKsiHeta -
		                     2 * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
		                     bsplinesHeta.BSPLineValues[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHeta, 2) +
		                     2 * bsplinesHeta.BSPLineValues[indexHeta, j] * Math.Pow(sumKsidHeta, 2) / Math.Pow(sumKsiHeta, 3));

	                    NurbsSecondDerivativeValueKsiHeta[k, i * supportHeta + j] =
		                    controlPoints[k].WeightFactor *
		                    (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
		                     bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] / sumKsiHeta -
		                     bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
		                     sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
		                     bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
		                     sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
		                     bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
		                     sumdKsidHeta / Math.Pow(sumKsiHeta, 2) +
		                     2 * bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
		                     sumdKsiHeta * sumKsidHeta / Math.Pow(sumKsiHeta, 3));
                    }
                }
            }
        }

		public NURBS2D(Element element, IList<ControlPoint> controlPoints, IVector parametricGaussPointKsi, IVector parametricGaussPointHeta)
		{
			var parametricPointsCount = parametricGaussPointKsi.Length * parametricGaussPointHeta.Length;

			BSPLines1D bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bsplinesHeta = new BSPLines1D(element.Patch.DegreeHeta, element.Patch.KnotValueVectorHeta, parametricGaussPointHeta);
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();


			int supportKsi = parametricGaussPointKsi.Length;
			int supportHeta = parametricGaussPointHeta.Length;
			int numberOfElementControlPoints = (element.Patch.DegreeKsi + 1) * (element.Patch.DegreeHeta + 1);

			NurbsValues = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);
			NurbsDerivativeValuesKsi = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);
			NurbsDerivativeValuesHeta = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);
			NurbsSecondDerivativeValueKsi = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);
			NurbsSecondDerivativeValueHeta = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);
			NurbsSecondDerivativeValueKsiHeta = new Matrix2D(numberOfElementControlPoints, parametricPointsCount);

			for (int i = 0; i < supportKsi; i++)
			{
				for (int j = 0; j < supportHeta; j++)
				{
					double sumKsiHeta = 0;
					double sumdKsiHeta = 0;
					double sumKsidHeta = 0;
					double sumdKsidKsi = 0;
					double sumdHetadHeta = 0;
					double sumdKsidHeta = 0;

					for (int k = 0; k < numberOfElementControlPoints; k++)
					{
						// Why type casting is needed.?

						int indexKsi = element.ControlPoints[k].ID / element.Patch.NumberOfControlPointsHeta;
						int indexHeta = element.ControlPoints[k].ID % element.Patch.NumberOfControlPointsHeta;
						sumKsiHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
									  bsplinesHeta.BSPLineValues[indexHeta, j] *
									  controlPoints[k].WeightFactor;
						sumdKsiHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
									   bsplinesHeta.BSPLineValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumKsidHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
									   bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumdKsidKsi += bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] *
									   bsplinesHeta.BSPLineValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumdHetadHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
										 bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] *
										 controlPoints[k].WeightFactor;
						sumdKsidHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
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

						NurbsSecondDerivativeValueKsi[k, i * supportHeta + j] =
							bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] / sumKsiHeta -
							 2 * bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * Math.Pow(sumdKsiHeta, 2) / Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueHeta[k, i * supportHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] / sumKsiHeta -
							 2 * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesHeta.BSPLineValues[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesHeta.BSPLineValues[indexHeta, j] * Math.Pow(sumKsidHeta, 2) / Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueKsiHeta[k, i * supportHeta + j] =
							controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
							 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] / sumKsiHeta -
							 bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
							 sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
							 sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
							 sumdKsidHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
							 sumdKsiHeta * sumKsidHeta / Math.Pow(sumKsiHeta, 3));
					}
				}
			}
		}

		public NURBS2D(Element element, IList<ControlPoint> controlPoints, Face face)
		{
			var degreeKsi = face.Degrees[0];
			var degreeHeta = face.Degrees[1];
			var knotValueVectorKsi = face.KnotValueVectors[0];
			var knotValueVectorHeta = face.KnotValueVectors[1];
			var numberOfControlPointsHeta = knotValueVectorHeta.Length - degreeHeta - 1;

			GaussQuadrature gauss = new GaussQuadrature();
			IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, element.Knots);

			IVector parametricGaussPointKsi = new Vector(degreeKsi + 1);
			for (int i = 0; i < degreeKsi + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i * (degreeHeta + 1)].Ksi;
			}

			IVector parametricGaussPointHeta = new Vector(degreeHeta + 1);
			for (int i = 0; i < degreeHeta + 1; i++)
			{
				parametricGaussPointHeta[i] = gaussPoints[i].Heta;
			}

			BSPLines1D bsplinesKsi = new BSPLines1D(degreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta, parametricGaussPointHeta);
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();


			int supportKsi = degreeKsi + 1;
			int supportHeta = degreeHeta + 1;
			int numberOfElementControlPoints = supportKsi * supportHeta;

			NurbsValues = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
			NurbsDerivativeValuesKsi = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
			NurbsDerivativeValuesHeta = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
			NurbsSecondDerivativeValueKsi = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
			NurbsSecondDerivativeValueHeta = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);
			NurbsSecondDerivativeValueKsiHeta = new Matrix2D(numberOfElementControlPoints, gaussPoints.Count);

			for (int i = 0; i < supportKsi; i++)
			{
				for (int j = 0; j < supportHeta; j++)
				{
					double sumKsiHeta = 0;
					double sumdKsiHeta = 0;
					double sumKsidHeta = 0;
					double sumdKsidKsi = 0;
					double sumdHetadHeta = 0;
					double sumdKsidHeta = 0;

					for (int k = 0; k < numberOfElementControlPoints; k++)
					{
						// Why type casting is needed.?

						int indexKsi = face.ControlPointsDictionary.First(cp=>cp.Value.ID== element.ControlPoints[k].ID).Key / numberOfControlPointsHeta;
						int indexHeta = face.ControlPointsDictionary.First(cp => cp.Value.ID == element.ControlPoints[k].ID).Key % numberOfControlPointsHeta;
						sumKsiHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
									  bsplinesHeta.BSPLineValues[indexHeta, j] *
									  controlPoints[k].WeightFactor;
						sumdKsiHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
									   bsplinesHeta.BSPLineValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumKsidHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
									   bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumdKsidKsi += bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] *
									   bsplinesHeta.BSPLineValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumdHetadHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
										 bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] *
										 controlPoints[k].WeightFactor;
						sumdKsidHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
										bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
										controlPoints[k].WeightFactor;
					}
					for (int k = 0; k < numberOfElementControlPoints; k++)
					{
						int indexKsi = face.ControlPointsDictionary.First(cp => cp.Value.ID == element.ControlPoints[k].ID).Key / numberOfControlPointsHeta;
						int indexHeta = face.ControlPointsDictionary.First(cp => cp.Value.ID == element.ControlPoints[k].ID).Key % numberOfControlPointsHeta;

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

						NurbsSecondDerivativeValueKsi[k, i * supportHeta + j] =
							bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] / sumKsiHeta -
							 2 * bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * Math.Pow(sumdKsiHeta, 2) / Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueHeta[k, i * supportHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] / sumKsiHeta -
							 2 * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesHeta.BSPLineValues[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesHeta.BSPLineValues[indexHeta, j] * Math.Pow(sumKsidHeta, 2) / Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueKsiHeta[k, i * supportHeta + j] =
							controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
							 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] / sumKsiHeta -
							 bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
							 sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
							 sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
							 sumdKsidHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
							 sumdKsiHeta * sumKsidHeta / Math.Pow(sumKsiHeta, 3));
					}
				}
			}
		}
	}
}
