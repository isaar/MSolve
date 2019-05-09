using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.Entities
{
	public class CollocationPoint: NaturalPoint, INode
    {
        private int _id;
        private bool _isBoundary;
        public int ID => _id;
        int INode.ID { get => _id; set => _id = value; }

        public bool IsBoundary
        {
            get => _isBoundary;
            set => _isBoundary = value;
        }

        public double X { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public double Y { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public double Z { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public List<Constraint> Constraints => new List<Constraint>();

        public Dictionary<int, ISubdomain> SubdomainsDictionary => throw new NotImplementedException();

        public CollocationPoint2D(int id, double xi, double eta, bool isBoundary=false) : base(xi, eta)
		{
            _id = id;
            _isBoundary = isBoundary;
        }

		public CollocationPoint2D(int id, double[] coordinates, bool isBoundary=false) : base(coordinates)
		{
            _id = id;
            _isBoundary = isBoundary;
        }
	}
}
