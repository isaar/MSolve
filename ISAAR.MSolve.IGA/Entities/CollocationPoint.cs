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
        public int ID => _id;
        int INode.ID { get => _id; }
        public double X { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public double Y { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public double Z { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public List<Constraint> Constraints => new List<Constraint>();

        public Dictionary<int, ISubdomain> SubdomainsDictionary => throw new NotImplementedException();

        public CollocationPoint(int id, double xi, double eta) : base(xi, eta)
		{
            _id = id;
		}

		public CollocationPoint(int id, double[] coordinates) : base(coordinates)
		{
            _id = id;
		}

        public int CompareTo(INode other) => this.ID - other.ID;
    }
}
