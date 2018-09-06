using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Problems.Structural.Elements;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
	public class Nurbs3D
	{
		private List<ControlPoint> ElementControlPoints()
		{
			return new List<ControlPoint>()
			{
				new ControlPoint(){ID=0,X =-1 ,Y=0.0,Z=0.0,WeightFactor =1.0 },
				new ControlPoint(){ID=1,X = -1,Y=0.0,Z=1.33333333333333,WeightFactor =1.0 },
				new ControlPoint(){ID=2,X =-1 ,Y=0.0,Z=2.33333333333333,WeightFactor =1.0 },
				new ControlPoint(){ID=3,X =-1,Y=0.0,Z=3.66666666666667,WeightFactor =1.0 },
				new ControlPoint(){ID=4,X =-2.0 ,Y=0.0,Z=0.0,WeightFactor =1.0 },
				new ControlPoint(){ID=5,X = -2.0,Y=0.0,Z=1.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=6,X = -2.0,Y=0.0,Z=2.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=7,X = -2.0,Y=0.0,Z=3.66666666666667,WeightFactor = 1.0},
				new ControlPoint(){ID=8,X = -3.0,Y=0.0,Z=0.0,WeightFactor =1.0 },
				new ControlPoint(){ID=9,X = -3.0,Y=0.0,Z=1.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=10,X = -3.0,Y=0.0,Z=2.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=11,X =-3 ,Y=0.0,Z=3.66666666666667,WeightFactor = 1.0},
				new ControlPoint(){ID=12,X = -4.0,Y=0.0,Z=0.0,WeightFactor = 1.0},
				new ControlPoint(){ID=13,X = -4.0,Y=0.0,Z=1.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=14,X = -4.0,Y=0.0,Z=2.33333333333333,WeightFactor =1.0 },
				new ControlPoint(){ID=15,X = -4,Y=0.0,Z=3.66666666666667,WeightFactor = 1.0},

				new ControlPoint(){ID=16,X = -0.902368926666667,Y=0.235702258881312,Z=0.0,WeightFactor = 0.902368926666667},
				new ControlPoint(){ID=17,X = -0.902368926666667,Y=0.235702258881312,Z=1.20315856888889,WeightFactor = 0.902368926666667},
				new ControlPoint(){ID=18,X = -0.902368926666667,Y=0.235702258881312,Z=2.10552749555556,WeightFactor = 0.902368926666667},
				new ControlPoint(){ID=19,X = -0.902368926666667,Y=0.235702258881312,Z=3.30868606444444,WeightFactor = 0.902368926666667},
				new ControlPoint(){ID=20,X = -1.93984208666667,Y=0.501985726293771,Z=0.0,WeightFactor = 0.967456308888889},
				new ControlPoint(){ID=21,X = -1.93984208666667,Y=0.501985726293771,Z=1.28994174518519,WeightFactor = 0.967456308888889},
				new ControlPoint(){ID=22,X = -1.93984208666667,Y=0.501985726293771,Z=2.25739805407407,WeightFactor = 0.967456308888889},
				new ControlPoint(){ID=23,X = -1.93984208666667,Y=0.501985726293771,Z=3.54733979925926,WeightFactor = 0.967456308888889},
				new ControlPoint(){ID=24,X = -2.97238577777778,Y=1.31230719555556,Z=0.0,WeightFactor = 1.0},
				new ControlPoint(){ID=25,X = -2.97238577777778,Y=1.31230719555556,Z=1.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=26,X = -2.97238577777778,Y=1.31230719555556,Z=2.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=27,X = -2.97238577777778,Y=1.31230719555556,Z=3.66666666666667,WeightFactor = 1.0},
				new ControlPoint(){ID=28,X = -4.0,Y=2.66666666666667,Z=0.0,WeightFactor =1.0 },
				new ControlPoint(){ID=29,X = -4.0,Y=2.66666666666667,Z=1.33333333333333,WeightFactor =1.0 },
				new ControlPoint(){ID=30,X = -4.0,Y=2.66666666666667,Z=2.33333333333333,WeightFactor =1.0},
				new ControlPoint(){ID=31,X = -4.0,Y=2.66666666666667,Z=3.66666666666667,WeightFactor =1.0},

				new ControlPoint(){ID=32,X = -0.770220056386995,Y=0.436886721934974,Z=0.0,WeightFactor = 0.85355339},
				new ControlPoint(){ID=33,X = -0.770220056386995,Y=0.436886721934974,Z=1.13807118666667,WeightFactor = 0.85355339 },
				new ControlPoint(){ID=34,X = -0.770220056386995,Y=0.436886721934974,Z=1.99162457666667,WeightFactor = 0.85355339 },
				new ControlPoint(){ID=35,X = -0.770220056386995,Y=0.436886721934974,Z=3.12969576333333,WeightFactor = 0.85355339},
				new ControlPoint(){ID=36,X = -1.71926689324011,Y=0.943474823978324,Z=0.0,WeightFactor = 0.951184463333333},
				new ControlPoint(){ID=37,X = -1.71926689324011,Y=0.943474823978324,Z=1.26824595111111,WeightFactor = 0.951184463333333},
				new ControlPoint(){ID=38,X = -1.71926689324011,Y=0.943474823978324,Z=2.21943041444444,WeightFactor = 0.951184463333333},
				new ControlPoint(){ID=39,X = -1.71926689324011,Y=0.943474823978324,Z=3.48767636555555,WeightFactor = 0.951184463333333},
				new ControlPoint(){ID=40,X = -2.79586020777778,Y=2.13117925,Z=0.0,WeightFactor = 1.0},
				new ControlPoint(){ID=41,X = -2.79586020777778,Y=2.13117925,Z=1.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=42,X = -2.79586020777778,Y=2.13117925,Z=2.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=43,X = -2.79586020777778,Y=2.13117925,Z=3.66666666666667,WeightFactor = 1.0},
				new ControlPoint(){ID=44,X = -4.0,Y=4.0,Z=0.0,WeightFactor =1.0 },
				new ControlPoint(){ID=45,X = -4.0,Y=4.0,Z=1.33333333333333,WeightFactor =1.0},
				new ControlPoint(){ID=46,X = -4.0,Y=4.0,Z=2.33333333333333,WeightFactor =1.0},
				new ControlPoint(){ID=47,X = -4.0,Y=4.0,Z=3.66666666666667,WeightFactor =1.0 },

				new ControlPoint(){ID=48,X = -0.436886721934974,Y=0.770220056386995,Z=0.0,WeightFactor = 0.85355339},
				new ControlPoint(){ID=49,X = -0.436886721934974,Y=0.770220056386995,Z=1.13807118666667,WeightFactor = 0.85355339},
				new ControlPoint(){ID=50,X = -0.436886721934974,Y=0.770220056386995,Z=1.99162457666667,WeightFactor = 0.85355339},
				new ControlPoint(){ID=51,X = -0.436886721934974,Y=0.770220056386995,Z=3.12969576333333,WeightFactor = 0.85355339},
				new ControlPoint(){ID=52,X = -0.957281946200547,Y=1.705459762129,Z=0.0,WeightFactor = 0.951184463333333},
				new ControlPoint(){ID=53,X = -0.957281946200547,Y=1.705459762129,Z=1.26824595111111,WeightFactor = 0.951184463333333},
				new ControlPoint(){ID=54,X = -0.957281946200547,Y=1.705459762129,Z=2.21943041444444,WeightFactor = 0.951184463333333},
				new ControlPoint(){ID=55,X = -0.957281946200547,Y=1.705459762129,Z=3.48767636555555,WeightFactor = 0.951184463333333},
				new ControlPoint(){ID=56,X = -2.14498637222222,Y=2.78205307666667,Z=0.0,WeightFactor =1.0 },
				new ControlPoint(){ID=57,X = -2.14498637222222,Y=2.78205307666667,Z=1.33333333333333,WeightFactor =1.0 },
				new ControlPoint(){ID=58,X = -2.14498637222222,Y=2.78205307666667,Z=2.33333333333333,WeightFactor =1.0 },
				new ControlPoint(){ID=59,X = -2.14498637222222,Y=2.78205307666667,Z=3.66666666666667,WeightFactor =1.0 },
				new ControlPoint(){ID=60,X = -4.0,Y=4.0,Z=0.0,WeightFactor = 1.0},
				new ControlPoint(){ID=61,X = -4.0,Y=4.0,Z=1.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=62,X = -4.0,Y=4.0,Z=2.33333333333333,WeightFactor = 1.0},
				new ControlPoint(){ID=63,X = -4.0,Y=4.0,Z=3.66666666666667,WeightFactor = 1.0}
			};
		}

		private List<Knot> ElementKnot()
		{
			return new List<Knot>()
			{
				new Knot() {ID = 0, Ksi = 0.0, Heta = 0.0, Zeta = 0.0},
				new Knot() {ID = 1, Ksi = 0.0, Heta = 0.0, Zeta = 0.166666666666667},
				new Knot() {ID = 2, Ksi = 0.0, Heta = 1.0, Zeta = 0.0},
				new Knot() {ID = 3, Ksi = 0.0, Heta = 1.0, Zeta = 0.166666666666667},
				new Knot() {ID = 4, Ksi = 0.5, Heta = 0.0, Zeta = 0.0},
				new Knot() {ID = 5, Ksi = 0.5, Heta = 0.0, Zeta = 0.166666666666667},
				new Knot() {ID = 6, Ksi = 0.5, Heta = 1.0, Zeta = 0.0},
				new Knot() {ID = 7, Ksi = 0.5, Heta = 1.0, Zeta = 0.166666666666667},
			};
		}

		private Vector KnotValueVectorKsi()
		{
			return new Vector(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0 });
		}

		private Vector KnotValueVectorHeta()
		{
			return new Vector(new double[]{ 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 });
		}

		private Vector KnotValueVectorZeta()
		{
			return new Vector(new double[]
			{
				0.0, 0.0, 0.0, 0.0, 0.166666666666667, 0.166666666666667, 0.333333333333333, 0.333333333333333, 0.5,
				0.5, 0.5, 0.666666666666667, 0.666666666666667, 0.833333333333333, 0.833333333333333, 1.0, 1.0, 1.0, 1.0
			});
		}

		[Fact]
		public void TestShapeNurbs3DPartitionOfUnity()
		{
			const double tolerance = 1e-10;
			var element = new NURBSElement3D();
			var patch = new Patch();
			foreach (var controlPoint in ElementControlPoints())
				element.ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
			foreach (var knot in ElementKnot())
				element.KnotsDictionary.Add(knot.ID, knot);
			patch.DegreeKsi = 3;
			patch.DegreeHeta = 3;
			patch.DegreeZeta = 3;
			patch.NumberOfControlPointsHeta = 4;
			patch.NumberOfControlPointsZeta = 15;
			patch.KnotValueVectorKsi = KnotValueVectorKsi();
			patch.KnotValueVectorHeta = KnotValueVectorHeta();
			patch.KnotValueVectorZeta = KnotValueVectorZeta();
			element.Patch = patch;

			var nurbs3D = new NURBS3D(element, element.ControlPoints);

			for (var p = 0; p < nurbs3D.NurbsValues.Columns; p++)
			{
				var sum = 0.0;
				for (var f = 0; f < nurbs3D.NurbsValues.Rows; f++)
					sum += nurbs3D.NurbsValues[f, p];
				Assert.True(Utilities.AreValuesEqual(1.0, sum, tolerance));
			}
		}

	}
}
