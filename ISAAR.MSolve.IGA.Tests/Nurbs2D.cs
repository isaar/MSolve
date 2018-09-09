using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Problems.Structural.Elements;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using Xunit;
using System.Linq;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using Moq;

namespace ISAAR.MSolve.IGA.Tests
{
	public class Nurbs2D
	{
		private List<ControlPoint> ElementControlPoints()
		{
			return new List<ControlPoint>
			{
				new ControlPoint {ID = 0, X = 0.0, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 1, X = 0.0, Y = 1.0/27.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 2, X = 0.0, Y = 1.0/9.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 3, X = 0.0, Y = 2.0/9.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 4, X = 1.0/27.0, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 5, X = 1.0/27.0, Y = 1.0/27.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 6, X = 1.0/27.0, Y = 1.0/9.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 7, X = 1.0/27.0, Y = 2.0/9.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 8, X = 1.0/9.0, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 9, X = 1.0/9.0, Y = 1.0/27.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 10, X = 1.0/9.0, Y = 1.0/9.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 11, X = 1.0/9.0, Y = 2.0/9.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 12, X = 2.0/9.0, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 13, X = 2.0/9.0, Y = 1.0/27.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 14, X = 2.0/9.0, Y = 1.0/9.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 15, X = 2.0/9.0, Y = 2.0/9.0, Z = 0.0, WeightFactor = 1.0},
			};
		}

		private List<Knot> ElementKnot()
		{
			return new List<Knot>()
			{
				new Knot(){ID=0,Ksi=0.0,Heta=0.0,Zeta =0.0 },
				new Knot(){ID=1,Ksi=0.0,Heta=1/9.0,Zeta =0.0 },
				new Knot(){ID=2,Ksi=1/9.0,Heta=0.0,Zeta =0.0 },
				new Knot(){ID=3,Ksi=1/9.0,Heta=1/9.0,Zeta =0.0 }
			};
		}

		private Vector KnotValueVector()
		{
			return new Vector(new double[16]
			{
				0.0, 0.0, 0.0, 0.0, 1/9.0, 2/9.0, 3/9.0, 4/9.0, 5/9.0, 6/9.0, 7/9.0,
				8/9.0, 1.0, 1.0, 1.0, 1.0
			});
		}

		[Fact]
		public void TestShapeNurbs2DPartitionOfUnity()
		{
			const double tolerance = 1e-10;
			var element = new NURBSElement2D();
			var patch = new Patch();
			foreach (var controlPoint in ElementControlPoints())
				element.ControlPointsDictionary.Add(controlPoint.ID,controlPoint);
			foreach (var knot in ElementKnot())
				element.KnotsDictionary.Add(knot.ID,knot);
			patch.DegreeKsi = 3;
			patch.DegreeHeta = 3;
			patch.NumberOfControlPointsHeta = 12;
			patch.KnotValueVectorKsi = KnotValueVector();
			patch.KnotValueVectorHeta = KnotValueVector();
			element.Patch = patch;



			var nurbs2D = new NURBS2D(element,element.ControlPoints);
			
			for (var p = 0; p < 16; p++)
			{
				var sum = 0.0;
				for (var f = 0; f < nurbs2D.NurbsValues.Rows; f++)
					sum += nurbs2D.NurbsValues[f, p];
				Assert.True(Utilities.AreValuesEqual(1.0, sum, tolerance));
			}

		}


	}
}
