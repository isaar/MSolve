using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface INode: IDiscretePoint
    {
		int ID { get; set; }
		double X { get; set; }
		double Y { get; set; }
		double Z { get; set; }

        List<Constraint> Constraints { get; }
        Dictionary<int, ISubdomain_v2> SubdomainsDictionary { get; }
    }

    public static class NodeExtensions
    {
        //TODO: perhaps this should be in the geometry project.
        public static double CalculateEuclidianDistanceFrom(this INode node1, INode node2)
        {
            double dx = node1.X - node2.X;
            double dy = node1.Y - node2.Y;
            double dz = node1.Z - node2.Z;
            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }
    }
}
