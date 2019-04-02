using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.Problems.SupportiveClasses
{
	public class NURBS3D
	{
		public double[,] NurbsValues { get; private set; }
		public double[,] NurbsDerivativeValuesKsi { get; private set; }
		public double[,] NurbsDerivativeValuesHeta { get; private set; }
		public double[,] NurbsDerivativeValuesZeta { get; private set; }
		public double[,] NurbsSecondDerivativeValueKsi { get; private set; }
		public double[,] NurbsSecondDerivativeValueHeta { get; private set; }
		public double[,] NurbsSecondDerivativeValueZeta { get; private set; }
		public double[,] NurbsSecondDerivativeValueKsiHeta { get; private set; }
		public double[,] NurbsSecondDerivativeValueKsiZeta { get; private set; }
		public double[,] NurbsSecondDerivativeValueHetaZeta { get; private set; }

		public NURBS3D(int numberOfControlPointsKsi, int numberOfControlPointsHeta, int numberOfControlPointsZeta,
			int degreeKsi, int degreeHeta, int degreeZeta, double[] knotValueVectorKsi,
			double[] knotValueVectorHeta, double[] knotValueVectorZeta, ControlPoint[] controlPoints,
			NaturalPoint3D collocationPoint)
		{
			BSPLines1D bsplinesKsi =
				new BSPLines1D(degreeKsi, knotValueVectorKsi, new double[] { collocationPoint.Xi });
			BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta,
				new double[] { collocationPoint.Eta });
			BSPLines1D bsplinesZeta = new BSPLines1D(degreeZeta, knotValueVectorZeta,
				new double[] { collocationPoint.Zeta });
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();
			bsplinesZeta.calculateBSPLinesAndDerivatives();

			int supportKsi = degreeKsi + 1;
			int supportHeta = degreeHeta + 1;
			int supportZeta = degreeZeta + 1;
			int numberOfGPKsi = 1;
			int numberOfGPHeta = 1;
			int numberOfGPZeta = 1;
			int numberOfElementControlPoints = supportKsi * supportHeta * supportZeta;

			NurbsValues = new double[numberOfElementControlPoints, 1];
			NurbsDerivativeValuesKsi = new double[numberOfElementControlPoints, 1];
			NurbsDerivativeValuesHeta = new double[numberOfElementControlPoints, 1];
			NurbsDerivativeValuesZeta = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueKsi = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueHeta = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueZeta = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueKsiHeta = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueKsiZeta = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueHetaZeta = new double[numberOfElementControlPoints, 1];

			for (int i = 0; i < numberOfGPKsi; i++)
			{
				for (int j = 0; j < numberOfGPHeta; j++)
				{
					for (int k = 0; k < numberOfGPZeta; k++)
					{
						double sumKsiHetaZeta = 0;
						double sumdKsiHetaZeta = 0;
						double sumKsidHetaZeta = 0;
						double sumKsiHetadZeta = 0;

						double sumdKsidKsi = 0;
						double sumdHetadHeta = 0;
						double sumdZetadZeta = 0;
						double sumdKsidHeta = 0;
						double sumdKsidZeta = 0;
						double sumdHetadZeta = 0;

						for (int m = 0; m < numberOfElementControlPoints; m++)
						{
							int indexKsi = controlPoints[m].ID /
										   (numberOfControlPointsHeta *
											numberOfControlPointsZeta);
							int indexHeta = controlPoints[m].ID %
											(numberOfControlPointsHeta *
											 numberOfControlPointsZeta) /
											numberOfControlPointsZeta;
							int indexZeta = controlPoints[m].ID %
											(numberOfControlPointsHeta *
											 numberOfControlPointsZeta) %
											numberOfControlPointsZeta;

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

							sumdKsidKsi += bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] *
							               bsplinesHeta.BSPLineValues[indexHeta, j] *
							               bsplinesZeta.BSPLineValues[indexZeta, k] *
							               controlPoints[m].WeightFactor;

							sumdHetadHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
							                 bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] *
							                 bsplinesZeta.BSPLineValues[indexZeta, k] *
							                 controlPoints[m].WeightFactor;

							sumdZetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
							                 bsplinesHeta.BSPLineValues[indexHeta, j] *
							                 bsplinesZeta.BSPLineSecondDerivativeValues[indexZeta, k]*
							                 controlPoints[m].WeightFactor;

							sumdKsidHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
							                bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
							                bsplinesZeta.BSPLineValues[indexZeta, k] *
							                controlPoints[m].WeightFactor;

							sumdKsidZeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
							                bsplinesHeta.BSPLineValues[indexHeta, j] *
							                bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
							                controlPoints[m].WeightFactor;

							sumdHetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
							                 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
							                 bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
							                 controlPoints[m].WeightFactor;
						}

						for (int m = 0; m < numberOfElementControlPoints; m++)
						{
							int indexKsi = controlPoints[m].ID /
										   (numberOfControlPointsHeta *
											numberOfControlPointsZeta);
							int indexHeta = controlPoints[m].ID %
											(numberOfControlPointsHeta *
											 numberOfControlPointsZeta) /
											numberOfControlPointsZeta;
							int indexZeta = controlPoints[m].ID %
											(numberOfControlPointsHeta *
											 numberOfControlPointsZeta) %
											numberOfControlPointsZeta;

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

							NurbsSecondDerivativeValueKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
								(bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] / sumKsiHetaZeta -
								 2 * bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumdKsiHetaZeta /
								 Math.Pow(sumKsiHetaZeta, 2) -
								 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHetaZeta, 2) +
								 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * Math.Pow(sumdKsiHetaZeta, 2) /
								 Math.Pow(sumKsiHetaZeta, 3)) *
								bsplinesHeta.BSPLineValues[indexHeta, j] *
								bsplinesZeta.BSPLineValues[indexZeta, k] *
								controlPoints[m].WeightFactor;

							NurbsSecondDerivativeValueHeta[m, i * supportHeta * supportZeta + j * supportZeta + k]=
								bsplinesKsi.BSPLineValues[indexKsi, i]*
								(bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] / sumKsiHetaZeta -
								 2 * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsidHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
								 bsplinesHeta.BSPLineValues[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHetaZeta, 2) +
								 2 * bsplinesHeta.BSPLineValues[indexHeta, j] * Math.Pow(sumKsidHetaZeta, 2) / Math.Pow(sumKsiHetaZeta, 3))*
								bsplinesZeta.BSPLineValues[indexZeta, k] *
								controlPoints[m].WeightFactor;

							NurbsSecondDerivativeValueZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
								bsplinesKsi.BSPLineValues[indexKsi, i] *
								bsplinesHeta.BSPLineValues[indexHeta, j] *
								(bsplinesZeta.BSPLineSecondDerivativeValues[indexZeta, k] / sumKsiHetaZeta -
								 2 * bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] * sumKsiHetadZeta /
								 Math.Pow(sumKsiHetaZeta, 2) -
								 bsplinesZeta.BSPLineValues[indexZeta, k] * sumdZetadZeta /
								 Math.Pow(sumKsiHetaZeta, 2) +
								 2 * bsplinesZeta.BSPLineValues[indexZeta, k] * Math.Pow(sumKsiHetadZeta, 2) /
								 Math.Pow(sumKsiHetaZeta, 3)) *
								controlPoints[m].WeightFactor;

							NurbsSecondDerivativeValueKsiHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
								(bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
								 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] / sumKsiHetaZeta -
								 bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
								 bsplinesHeta.BSPLineValues[indexHeta, j] *
								 sumKsidHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
								 bsplinesKsi.BSPLineValues[indexKsi, i] *
								 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
								 sumdKsiHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
								 bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
								 sumdKsidHeta / Math.Pow(sumKsiHetaZeta, 2) +
								 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
								 sumdKsiHetaZeta * sumKsidHetaZeta / Math.Pow(sumKsiHetaZeta, 3)) *
								bsplinesZeta.BSPLineValues[indexZeta, k] *
								controlPoints[m].WeightFactor;

							NurbsSecondDerivativeValueKsiZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
								(bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
								 bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] / sumKsiHetaZeta -
								 bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
								 bsplinesZeta.BSPLineValues[indexZeta, k] *
								 sumKsiHetadZeta / Math.Pow(sumKsiHetaZeta, 2) -
								 bsplinesKsi.BSPLineValues[indexKsi, i] *
								 bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
								 sumdKsiHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
								 bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesZeta.BSPLineValues[indexZeta, k] *
								 sumdKsidZeta / Math.Pow(sumKsiHetaZeta, 2) +
								 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesZeta.BSPLineValues[indexZeta, k] *
								 sumdKsiHetaZeta * sumKsiHetadZeta / Math.Pow(sumKsiHetaZeta, 3)) *
								bsplinesHeta.BSPLineValues[indexHeta, j] *
								controlPoints[m].WeightFactor;

							NurbsSecondDerivativeValueHetaZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
								bsplinesKsi.BSPLineValues[indexKsi, i] *
								(bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
								 bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] / sumKsiHetaZeta -
								 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
								 bsplinesZeta.BSPLineValues[indexZeta, k] *
								 sumKsiHetadZeta / Math.Pow(sumKsiHetaZeta, 2) -
								 bsplinesHeta.BSPLineValues[indexHeta, j] *
								 bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
								 sumKsidHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
								 bsplinesHeta.BSPLineValues[indexHeta, j] * bsplinesZeta.BSPLineValues[indexZeta, k] *
								 sumdHetadZeta / Math.Pow(sumKsiHetaZeta, 2) +
								 2 * bsplinesHeta.BSPLineValues[indexHeta, j] *
								 bsplinesZeta.BSPLineValues[indexZeta, k] *
								 sumKsidHetaZeta * sumKsiHetadZeta / Math.Pow(sumKsiHetaZeta, 3)) *
								controlPoints[m].WeightFactor;


						}
					}
				}
			}
		}


		public NURBS3D(int numberOfControlPointsKsi, int numberOfControlPointsHeta, int numberOfControlPointsZeta,
			int degreeKsi, int degreeHeta, int degreeZeta, double[] knotValueVectorKsi,
			double[] knotValueVectorHeta, double[] knotValueVectorZeta, ControlPoint[] controlPoints,
			GaussLegendrePoint3D[] gaussPoints)
		{
			var numberOfGaussPoints = gaussPoints.Length;
			var parametricGaussPointKsi = new double[degreeKsi + 1];
			for (int i = 0; i < degreeKsi + 1; i++)
				parametricGaussPointKsi[i] = gaussPoints[i * (degreeZeta + 1) * (degreeHeta + 1)].Ksi;

			var parametricGaussPointHeta = new double[degreeHeta + 1];
			for (int i = 0; i < degreeHeta + 1; i++)
				parametricGaussPointHeta[i] = gaussPoints[i * (degreeZeta + 1)].Heta;

			var parametricGaussPointZeta = new double[degreeZeta + 1];
			for (int i = 0; i < degreeZeta + 1; i++)
				parametricGaussPointZeta[i] = gaussPoints[i].Zeta;

			BSPLines1D bsplinesKsi =
				new BSPLines1D(degreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta,
				parametricGaussPointHeta);
			BSPLines1D bsplinesZeta = new BSPLines1D(degreeZeta, knotValueVectorZeta,
				parametricGaussPointZeta);
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();
			bsplinesZeta.calculateBSPLinesAndDerivatives();

			int supportKsi = degreeKsi + 1;
			int supportHeta = degreeHeta + 1;
			int supportZeta = degreeZeta + 1;
			int numberOfElementControlPoints = supportKsi * supportHeta * supportZeta;

			NurbsValues = new double[numberOfElementControlPoints, numberOfGaussPoints];
			NurbsDerivativeValuesKsi = new double[numberOfElementControlPoints, numberOfGaussPoints];
			NurbsDerivativeValuesHeta = new double[numberOfElementControlPoints, numberOfGaussPoints];
			NurbsDerivativeValuesZeta = new double[numberOfElementControlPoints, numberOfGaussPoints];

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
							int indexKsi = controlPoints[m].ID /
							               (numberOfControlPointsHeta *
							                numberOfControlPointsZeta);
							int indexHeta = controlPoints[m].ID %
							                (numberOfControlPointsHeta *
							                 numberOfControlPointsZeta) /
							                numberOfControlPointsZeta;
							int indexZeta = controlPoints[m].ID %
							                (numberOfControlPointsHeta *
							                 numberOfControlPointsZeta) %
							                numberOfControlPointsZeta;

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
							int indexKsi = controlPoints[m].ID /
							               (numberOfControlPointsHeta *
							                numberOfControlPointsZeta);
							int indexHeta = controlPoints[m].ID %
							                (numberOfControlPointsHeta *
							                 numberOfControlPointsZeta) /
							                numberOfControlPointsZeta;
							int indexZeta = controlPoints[m].ID %
							                (numberOfControlPointsHeta *
							                 numberOfControlPointsZeta) %
							                numberOfControlPointsZeta;

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

		public NURBS3D(int numberOfControlPointsKsi, int numberOfControlPointsHeta, int numberOfControlPointsZeta,
			int degreeKsi, int degreeHeta, int degreeZeta, double[] knotValueVectorKsi,
			double[] knotValueVectorHeta, double[] knotValueVectorZeta, ControlPoint[] controlPoints, double[] parametricGaussPointKsi,
			double[] parametricGaussPointHeta, double[] parametricGaussPointZeta)
		{
			var parametricPointsCount = parametricGaussPointKsi.Length * parametricGaussPointHeta.Length *
			                            parametricGaussPointZeta.Length;

			BSPLines1D bsplinesKsi = new BSPLines1D(degreeKsi, knotValueVectorKsi,
				parametricGaussPointKsi);
			BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta,
				parametricGaussPointHeta);
			BSPLines1D bsplinesZeta = new BSPLines1D(degreeZeta, knotValueVectorZeta,
				parametricGaussPointZeta);
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();
			bsplinesZeta.calculateBSPLinesAndDerivatives();

			int supportKsi = parametricGaussPointKsi.Length;
			int supportHeta = parametricGaussPointHeta.Length;
			int supportZeta = parametricGaussPointZeta.Length;
			int numberOfElementControlPoints = (degreeKsi + 1) * (degreeHeta + 1) *
			                                   (degreeZeta + 1);

			NurbsValues = new double[numberOfElementControlPoints, parametricPointsCount];
			NurbsDerivativeValuesKsi = new double[numberOfElementControlPoints, parametricPointsCount];
			NurbsDerivativeValuesHeta = new double[numberOfElementControlPoints, parametricPointsCount];
			NurbsDerivativeValuesZeta = new double[numberOfElementControlPoints, parametricPointsCount];

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
							int indexKsi = controlPoints[m].ID /
							               (numberOfControlPointsHeta *
							                numberOfControlPointsZeta);
							int indexHeta = controlPoints[m].ID %
							                (numberOfControlPointsHeta *
							                 numberOfControlPointsZeta) /
							                numberOfControlPointsZeta;
							int indexZeta = controlPoints[m].ID %
							                (numberOfControlPointsHeta *
							                 numberOfControlPointsZeta) %
							                numberOfControlPointsZeta;

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
							int indexKsi = controlPoints[m].ID /
							               (numberOfControlPointsHeta *
							                numberOfControlPointsZeta);
							int indexHeta = controlPoints[m].ID %
							                (numberOfControlPointsHeta *
							                 numberOfControlPointsZeta) /
							                numberOfControlPointsZeta;
							int indexZeta = controlPoints[m].ID %
							                (numberOfControlPointsHeta *
							                 numberOfControlPointsZeta) %
							                numberOfControlPointsZeta;

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