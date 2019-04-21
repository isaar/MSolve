using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.Entities
{
    public class CollocationPoint3D:NaturalPoint3D, INode
    {
        public CollocationPoint3D(int id,double xi, double eta, double zeta) : base(xi, eta, zeta)
        {
            _id = id;
        }

        public CollocationPoint3D(int id,double[] coordinates) : base(coordinates)
        {
            _id = id;
        }

        private int _id;
        public int ID => _id;

        int INode.ID { get => _id; set => _id = value; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public List<Constraint> Constraints { get; }
        public Dictionary<int, ISubdomain> SubdomainsDictionary { get; }
    }
}
