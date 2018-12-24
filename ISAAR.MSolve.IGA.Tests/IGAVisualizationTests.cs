using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.Numerical.LinearAlgebra;
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
			var spanKsi=ParaviewResults2D.FindSpan(numberOfControlPointsKsi, degreeKsi, coordinateKsi, knotValueVectorKsi);
			Assert.Equal(4, spanKsi);
		}

		[Fact]
		public void BasisFunctionsTest()
		{
			var spanKsi = ParaviewResults2D.FindSpan(numberOfControlPointsKsi, degreeKsi, coordinateKsi, knotValueVectorKsi);
			var basisFunctionsKsi = ParaviewResults2D.BasisFunctions(spanKsi, coordinateKsi, degreeKsi, knotValueVectorKsi);

			var expectedBasisFunctions=new Vector(new double[]{0.5,0.5,0});
			for (int i = 0; i < basisFunctionsKsi.Length; i++)
			{
				Utilities.AreValuesEqual(expectedBasisFunctions[i], basisFunctionsKsi[i], 1e-14);
			}

		}

		[Fact]
		public void SurfacePoint2DTest()
		{
			var point=ParaviewResults2D.SurfacePoint2D(numberOfControlPointsKsi, degreeKsi, knotValueVectorKsi,
				numberOfControlPointsHeta, degreeHeta, knotValueVectorHeta, projectiveCoordinates, coordinateKsi, coordinateHeta);

			var expectedPoint = new double[] {10, 0,0, 1};
			for (int i = 0; i < point.Length; i++)
			{
				Utilities.AreValuesEqual(expectedPoint[i], point[i], 1e-14);
			}
		}

	}
}
