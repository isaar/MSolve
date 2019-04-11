using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.Entities
{
	public class CollocationPoint2D:NaturalPoint2D, INode
    {
        private int _id;
        public int ID => _id;
        int INode.ID { get => _id; set => _id = value; }
        public double X { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public double Y { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public double Z { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public List<Constraint> Constraints => new List<Constraint>();

        public CollocationPoint2D(int id, double xi, double eta) : base(xi, eta)
		{
            _id = id;
		}

		public CollocationPoint2D(int id, double[] coordinates) : base(coordinates)
		{
            _id = id;
		}
	}
}
