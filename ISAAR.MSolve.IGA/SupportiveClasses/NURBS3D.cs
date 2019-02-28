using ISAAR.MSolve.IGA.Entities;
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
		public double[,] NurbsValues { get; private set; }
		public double[,] NurbsDerivativeValuesKsi { get; private set; }
		public double[,] NurbsDerivativeValuesHeta { get; private set; }
		public double[,] NurbsDerivativeValuesZeta { get; private set; }

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