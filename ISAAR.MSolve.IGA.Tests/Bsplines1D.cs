using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
	public class Bsplines1D
	{
		[Fact]
		public void TestBSplines1DValues()
		{
			const double tolerance = 1e-7;
			var degree = 3;
			var knotValueVector = new Vector(new double[]
			{
				0, 0, 0, 0, 0.11111111111111111, 0.22222222222222222, 0.33333333333333333, 0.44444444444444444,
				0.55555555555555555, 0.66666666666666666, 0.7777777777777777, 0.8888888888888888, 1, 1, 1, 1
			});

			var parametricCoordinates = new Vector(new double[]
				{0.007714649348171322, 0.0366677197641736, 0.0744433912358264, 0.10339646165182867});

			var bsplines1D= new BSPLines1D(degree, knotValueVector,parametricCoordinates);
			bsplines1D.calculateBSPLinesAndDerivatives();

			var expectedValues = new double[4, 4]
			{
				{0.805832094644764, 0.300750235878434, 0.0359400966193526, 0.000334715714594480},
				{0.187187770657043, 0.562845453025830, 0.516291631351937, 0.305103716393759},
				{0.00692434874576065, 0.130414294992511, 0.397643232715638, 0.560256218784186},
				{5.57859524324130e-05, 0.00599001610322543, 0.0501250393130722, 0.134305349107461},
			};

			for (var i = 0; i < 4; i++)
			{
				for (var j = 0; j < 4; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedValues[i, j], bsplines1D.BSPLineValues[i, j],
						tolerance));
				}
			}
		}

		[Fact]
		public void TestBSplines1DDerivativeValues()
		{
			const double tolerance = 1e-7;
			var degree = 3;
			var knotValueVector = new Vector(new double[]
			{
				0, 0, 0, 0, 0.11111111111111111, 0.22222222222222222, 0.33333333333333333, 0.44444444444444444,
				0.55555555555555555, 0.66666666666666666, 0.7777777777777777, 0.8888888888888888, 1, 1, 1, 1
			});

			var parametricCoordinates = new Vector(new double[]
				{0.007714649348171322, 0.0366677197641736, 0.0744433912358264, 0.10339646165182867});

			var bsplines1D = new BSPLines1D(degree, knotValueVector, parametricCoordinates);
			bsplines1D.calculateBSPLinesAndDerivatives();

			var expectedDerivativeValues = new double[4, 4]
			{
				{-23.3808414997539, -12.1199570808756, -2.94046890408451, -0.130161086714501},
				{21.6038025213095, 5.41505284733459, -6.05930737365432, -7.45954799498977},
				{1.75534546399199, 6.21482608286030, 6.97978343092623, 3.69290216507861},
				{0.0216935144524168, 0.490078150680752, 2.01999284681261, 3.89680691662566},
			};

			for (var i = 0; i < 4; i++)
			{
				for (var j = 0; j < 4; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedDerivativeValues[i, j], bsplines1D.BSPLineDerivativeValues[i, j],
						tolerance));
				}
			}
		}


		[Fact]
		public void TestBSplines1DPartitionOfUnity()
		{
			const double tolerance = 1e-10;
			var degree = 3;
			var knotValueVector = new Vector(new double[]
			{
				0, 0, 0, 0, 0.11111111111111111, 0.22222222222222222, 0.33333333333333333, 0.44444444444444444,
				0.55555555555555555, 0.66666666666666666, 0.7777777777777777, 0.8888888888888888, 1, 1, 1, 1
			});

			var parametricCoordinates = new Vector(new double[]
				{0.007714649348171322, 0.0366677197641736, 0.0744433912358264, 0.10339646165182867});

			var bsplines1D = new BSPLines1D(degree, knotValueVector, parametricCoordinates);
			bsplines1D.calculateBSPLinesAndDerivatives();
			
			for (var p = 0; p < 4; p++)
			{
				var sum = 0.0;
				for (var f = 0; f < bsplines1D.BSPLineValues.GetLength(0); f++)
					sum += bsplines1D.BSPLineValues[f, p];
				Assert.True(Utilities.AreValuesEqual(1.0, sum, tolerance));
			}
		}
	}
}
