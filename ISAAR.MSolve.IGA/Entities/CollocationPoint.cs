using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.Entities
{
	public class CollocationPoint2D:NaturalPoint2D
	{
		public int ID { get; }

		public CollocationPoint2D(int id, double xi, double eta) : base(xi, eta)
		{
			ID = id;
		}

		public CollocationPoint2D(int id, double[] coordinates) : base(coordinates)
		{
			ID = id;
		}
	}
}
