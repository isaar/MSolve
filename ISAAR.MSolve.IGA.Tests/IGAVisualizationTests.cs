using System;
using System.Collections.Generic;
using System.Diagnostics.Tracing;
using System.Text;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
	public class IGAVisualizationTests
	{
		private readonly int numberOfControlPointsKsi = 9;
		private readonly int degreeKsi = 2;
		private readonly double coordinateKsi = 0.25;
		private readonly Vector knotValueVectorKsi = new Vector(new double[]
			{0, 0, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1, 1});

		private readonly int numberOfControlPointsHeta = 5;
		private readonly int degreeHeta = 2;
		private readonly double coordinateHeta = 0.0;
		private readonly Vector knotValueVectorHeta = new Vector(new double[]
			{0,  0 ,  0,   0.25,   0.5,   0.75,   1,   1,   1});

		private readonly double[,] projectiveCoordinates = new double[,]
		{
			{0, 0,0, 1},
			{0, 1.25000000000000,0, 1},
			{0, 3.75000000000000,0, 1},
			{0, 6.25000000000000,0, 1},
			{0, 8.75000000000000, 0,1},
			{0, 10,0,1},

			{2.50000000000000, 0,0, 1},
			{2.50000000000000, 1.25000000000000,0,1},
			{2.50000000000000, 3.75000000000000,0, 1},
			{2.50000000000000, 6.25000000000000,0, 1},
			{2.50000000000000, 8.75000000000000,0, 1},
			{2.50000000000000, 10,0, 1},

			{7.50000000000000, 0,0, 1},
			{7.50000000000000, 1.25000000000000,0, 1},
			{7.50000000000000, 3.75000000000000,0, 1},
			{7.50000000000000, 6.25000000000000,0, 1},
			{7.50000000000000, 8.75000000000000,0, 1},
			{7.50000000000000, 10,0, 1},

			{12.5000000000000, 0,0, 1},
			{12.5000000000000, 1.25000000000000,0, 1},
			{12.5000000000000, 3.75000000000000,0, 1},
			{12.5000000000000, 6.25000000000000,0, 1},
			{12.5000000000000, 8.75000000000000,0, 1},
			{12.5000000000000, 10,0, 1},

			{17.5000000000000, 0,0, 1},
			{17.5000000000000, 1.25000000000000,0, 1},
			{17.5000000000000, 3.75000000000000,0, 1},
			{17.5000000000000, 6.25000000000000,0, 1},
			{17.5000000000000, 8.75000000000000,0, 1},
			{17.5000000000000, 10,0, 1},

			{22.5000000000000, 0,0, 1},
			{22.5000000000000, 1.25000000000000,0, 1},
			{22.5000000000000, 3.75000000000000,0, 1},
			{22.5000000000000, 6.25000000000000,0, 1},
			{22.5000000000000, 8.75000000000000,0, 1},
			{22.5000000000000, 10,0, 1},

			{27.5000000000000, 0,0, 1},
			{27.5000000000000, 1.25000000000000,0, 1},
			{27.5000000000000, 3.75000000000000,0, 1},
			{27.5000000000000, 6.25000000000000,0, 1},
			{27.5000000000000, 8.75000000000000,0, 1},
			{27.5000000000000, 10,0, 1},

			{32.5000000000000, 0,0, 1},
			{32.5000000000000, 1.25000000000000,0, 1},
			{32.5000000000000, 3.75000000000000,0, 1},
			{32.5000000000000, 6.25000000000000,0, 1},
			{32.5000000000000, 8.75000000000000,0, 1},
			{32.5000000000000, 10,0, 1},

			{37.5000000000000, 0, 0,1},
			{37.5000000000000, 1.25000000000000,0, 1},
			{37.5000000000000, 3.75000000000000,0, 1},
			{37.5000000000000, 6.25000000000000,0, 1},
			{37.5000000000000, 8.75000000000000,0, 1},
			{37.5000000000000, 10,0, 1},

			{40, 0,0, 1},
			{40, 1.25000000000000,0, 1},
			{40, 3.75000000000000,0, 1},
			{40, 6.25000000000000,0, 1},
			{40, 8.75000000000000,0, 1},
			{40, 10,0, 1}
		};


		[Fact]
		public void FindSpanTest()
		{
			var spanKsi=ParaviewNurbs2D.FindSpan(numberOfControlPointsKsi, degreeKsi, coordinateKsi, knotValueVectorKsi);
			Assert.Equal(4, spanKsi);
		}

		[Fact]
		public void BasisFunctionsTest()
		{
			var spanKsi = ParaviewNurbs2D.FindSpan(numberOfControlPointsKsi, degreeKsi, coordinateKsi, knotValueVectorKsi);
			var basisFunctionsKsi = ParaviewNurbs2D.BasisFunctions(spanKsi, coordinateKsi, degreeKsi, knotValueVectorKsi);

			var expectedBasisFunctions=new Vector(new double[]{0.5,0.5,0});
			for (int i = 0; i < basisFunctionsKsi.Length; i++)
			{
				Utilities.AreValuesEqual(expectedBasisFunctions[i], basisFunctionsKsi[i], 1e-14);
			}

		}

		[Fact]
		public void SurfacePoint2DTest()
		{
			var point=ParaviewNurbs2D.SurfacePoint2D(numberOfControlPointsKsi, degreeKsi, knotValueVectorKsi,
				numberOfControlPointsHeta, degreeHeta, knotValueVectorHeta, projectiveCoordinates, coordinateKsi, coordinateHeta);

			var expectedPoint = new double[] {10, 0,0, 1};
			for (int i = 0; i < point.Length; i++)
			{
				Utilities.AreValuesEqual(expectedPoint[i], point[i], 1e-14);
			}
		}

		private int numberOfCPKsi3D = 20;
		private int degreeKsi3D = 2;
		private Vector knotValueVectorKsi3D = new Vector(new double[]
		{
			0,0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,19,19
		});
		private double coordinateKsi3D = 0.052631578947368;

		private int numberOfCPHeta3D = 6;
		private int degreeHeta3D = 2;
		private Vector knotValueVectorHeta3D = new Vector(new double[]
		{
			0,0,0,1,2,3,4,5,5,5
		});
		private double coordinateHeta3D = 0.0;

		private int numberofCPZeta3D = 4;
		private int degreeZeta3D = 2;
		private Vector knotValueVectorZeta3D = new Vector(new double[]
		{
			0, 0, 0, 1, 2, 3, 3, 3
		});
		private double coordinateZeta3D = 0.0;

		[Fact]
		public void SolidPoint3D()
		{
			Matrix<double> projcoord = MatlabReader.Read<double>("..\\..\\..\\InputFiles\\SolidPointProjectiveCP.mat", "projcoord");

			var projectiveCPCoordinates = new double[735, 4];
			var counter = 0;
			for (int i = 0; i < numberOfCPKsi3D+1; i++)
			{
				for (int j = 0; j < numberOfCPHeta3D+1; j++)
				{
					for (int k = 0; k < numberofCPZeta3D+1; k++)
					{
						var index = k * (numberOfCPHeta3D + 1) * (numberOfCPKsi3D + 1) + j * (numberOfCPKsi3D + 1) + i;
						projectiveCPCoordinates[counter, 0] = projcoord.At(index, 0);
						projectiveCPCoordinates[counter, 1] = projcoord.At(index, 1);
						projectiveCPCoordinates[counter, 2] = projcoord.At(index, 2);
						projectiveCPCoordinates[counter, 3] = projcoord.At(index, 3);
					}
				}
			}


			var point = ParaviewNurbs3D.SolidPoint3D(numberOfCPKsi3D, degreeKsi3D, knotValueVectorKsi3D,
				numberOfCPHeta3D, degreeHeta3D, knotValueVectorHeta3D, numberofCPZeta3D, degreeZeta3D,
				knotValueVectorZeta3D, projectiveCPCoordinates, coordinateKsi3D, coordinateHeta3D, coordinateZeta3D);

			var expectedPoint = new double[] { 0.75, 0, 0, 1 };
			for (int i = 0; i < point.Length; i++)
			{
				Utilities.AreValuesEqual(expectedPoint[i], point[i], 1e-14);
			}
		}

	}
}
