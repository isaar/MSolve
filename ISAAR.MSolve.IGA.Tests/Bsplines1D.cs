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
		private const double Tolerance = 1e-14;

		private static Vector KnotValueVector
		{
			get
			{
				var knotValueVector = new Vector(new double[]
				{
					0, 0, 0, 0, 0.11111111111111111, 0.22222222222222222, 0.33333333333333333, 0.44444444444444444,
					0.55555555555555555, 0.66666666666666666, 0.7777777777777777, 0.8888888888888888, 1, 1, 1, 1
				});
				return knotValueVector;
			}
		}

		private static Vector ParametricCoordinates
		{
			get
			{
				return new Vector(new double[]
					{0.007714649348171322, 0.0366677197641736, 0.0744433912358264, 0.10339646165182867});
			}
		}


		[Fact]
		public void TestBSplines1DValues()
		{
			
			var degree = 3;
			var knotValueVector = KnotValueVector;

			var parametricCoordinates = ParametricCoordinates;

			var bsplines1D= new BSPLines1D(degree, knotValueVector,parametricCoordinates);
			bsplines1D.calculateBSPLinesAndDerivatives();

			var expectedValues = new double[4, 4]
			{
				{0.8058320948251374, 0.3007502363228445, 0.035940096838251105, 3.3471572805268126E-4},
				{0.18718777049037877, 0.5628454528272727, 0.5162916318030126, 0.30510371716505036},
				{0.006924348732218876, 0.13041429476462738, 0.39764323219603925, 0.5602562184023526},
				{5.5785952265056314E-5, 0.005990016085255392, 0.05012503916269711, 0.13430534870454433},
			};

			for (var i = 0; i < 4; i++)
			{
				for (var j = 0; j < 4; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedValues[i, j], bsplines1D.BSPLineValues[i, j],
						Tolerance));
				}
			}
		}
		
		[Fact]
		public void TestBSplines1DDerivativeValues()
		{
			var degree = 3;
			var knotValueVector = KnotValueVector;

			var parametricCoordinates = ParametricCoordinates;

			var bsplines1D = new BSPLines1D(degree, knotValueVector, parametricCoordinates);
			bsplines1D.calculateBSPLinesAndDerivatives();

			var expectedDerivativeValues = new double[4, 4]
			{
				{-23.380841503242923, -12.119957092815207, -2.940468916024088, -0.13016109020350003},
				{21.603802526477924, 5.4150528637737025, -6.0593073618049385, -7.4595480014466915},
				{1.7553454623559663, 6.214826079340906, 6.979783435056406, 3.6929021828181523},
				{0.02169351440903006, 0.49007814970059616, 2.01999284277262, 3.8968069088320396},
			};

			for (var i = 0; i < 4; i++)
			{
				for (var j = 0; j < 4; j++)
				{
					Assert.True(Utilities.AreValuesEqual(expectedDerivativeValues[i, j], bsplines1D.BSPLineDerivativeValues[i, j],
						Tolerance));
				}
			}
		}


		[Fact]
		public void TestBSplines1DPartitionOfUnity()
		{
			var degree = 3;
			var knotValueVector = KnotValueVector;

			var parametricCoordinates = ParametricCoordinates;

			var bsplines1D = new BSPLines1D(degree, knotValueVector, parametricCoordinates);
			bsplines1D.calculateBSPLinesAndDerivatives();
			
			for (var p = 0; p < 4; p++)
			{
				var sum = 0.0;
				for (var f = 0; f < bsplines1D.BSPLineValues.GetLength(0); f++)
					sum += bsplines1D.BSPLineValues[f, p];
				Assert.True(Utilities.AreValuesEqual(1.0, sum, Tolerance));
			}
		}
	}
}
