using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using Xunit;

namespace ISAAR.MSolve.FEM.Tests.Interpolation
{
	/// <summary>
	/// Unit testing implementations of <see cref="IIsoparametricInterpolation3D"/>
	/// </summary>
    public class IsoparametricInterpolation3D
	{
		private const int numRandomPoints = 10;
		private delegate NaturalPoint3D[] GenerateRandomPoints();

		public static readonly IEnumerable<object[]> interpolations = new List<object[]>()
		{
			new object[]{InterpolationTet4.UniqueInstance},
			new object[]{InterpolationTet10.UniqueInstance},
			new object[]{InterpolationHexa8.UniqueInstance}

		};

		private static readonly Dictionary<IIsoparametricInterpolation3D, GenerateRandomPoints> pointGenerators= 
			new Dictionary<IIsoparametricInterpolation3D, GenerateRandomPoints>()
			{
				{InterpolationTet4.UniqueInstance,GenerateRandomPointsInTetrahedron},
				{InterpolationTet10.UniqueInstance, GenerateRandomPointsInTetrahedron},
				{InterpolationHexa8.UniqueInstance, GenerateRandomPointsInCube }
			};


		[Theory]
		[MemberData(nameof(interpolations))]
		private static void TestPartitionOfUnity(IIsoparametricInterpolation3D interpolation)
		{
			double tolerance = 1e-10;
			NaturalPoint3D[] points = pointGenerators[interpolation]();
			for (int p = 0; p < points.Length; p++)
			{
				EvalShapeFunctions3D shapeFunctions = interpolation.EvaluateFunctionsAt(points[p]);
				double sum = 0.0;
				for (int f = 0; f < interpolation.NumFunctions; f++) sum += shapeFunctions[f];
				Assert.True(Utilities.AreValuesEqual(1.0,sum,tolerance));
			}
		}

		[Theory]
		[MemberData(nameof(interpolations))]
		private static void TestValuesAtNodes(IIsoparametricInterpolation3D interpolation)
		{
			double tolerance = 1e-10;
			for (int n = 0; n < interpolation.NodalNaturalCoordinates.Count; n++)
			{
				EvalShapeFunctions3D shapeFunctions = interpolation.EvaluateFunctionsAt(interpolation.NodalNaturalCoordinates[n]);
				for (int f = 0; f < interpolation.NumFunctions; f++)
				{
					if (f==n) Assert.True(Utilities.AreValuesEqual(1.0,shapeFunctions[f], tolerance));
					else Assert.True(Utilities.AreValuesEqual(0.0,shapeFunctions[f],tolerance));
				}
			}
		}


		/// <summary>
		/// Generates random point in the tetrahedron.
		/// </summary>
		/// <returns></returns>
		private static NaturalPoint3D[] GenerateRandomPointsInTetrahedron()
		{
			var rand= new Random();
			var randomPoints= new NaturalPoint3D[numRandomPoints];
			for (int i = 0; i < numRandomPoints; i++)
			{
				double xi = rand.NextDouble();
				double eta = rand.NextDouble() * xi;
				double zeta = rand.NextDouble() * eta;
				randomPoints[i]= new NaturalPoint3D(xi,eta,zeta);
			}

			return randomPoints;
		}

		/// <summary>
		/// Generates ranom points in the parent cube.
		/// </summary>
		/// <returns></returns>
		private static NaturalPoint3D[] GenerateRandomPointsInCube()
		{
			var rand = new Random();
			var randomPoints= new NaturalPoint3D[numRandomPoints];
			for (int i = 0; i < numRandomPoints; i++)
			{
				double xi = -1 + rand.NextDouble() * 2.0;
				double eta = -1 + rand.NextDouble() * 2.0;
				double zeta = -1 + rand.NextDouble() * 2.0;
				randomPoints[i] = new NaturalPoint3D(xi, eta, zeta);
			}

			return randomPoints;
		}
	}
}
