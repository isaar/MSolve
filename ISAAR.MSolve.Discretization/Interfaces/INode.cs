using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface INode : IComparable<INode>
    {
		int ID { get; set; }
		double X { get; set; }
		double Y { get; set; }
		double Z { get; set; }

        List<Constraint> Constraints { get; }
        Dictionary<int, ISubdomain> SubdomainsDictionary { get; }
    }
}
