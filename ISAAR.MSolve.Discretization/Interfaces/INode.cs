using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface INode
    {
		int ID { get; set; }
		double X { get; set; }
		double Y { get; set; }
		double Z { get; set; }

		List<DOFType> Constraints { get; }
		Dictionary<int, IElement> ElementsDictionary { get; }
		Dictionary<int, ISubdomain> SubdomainsDictionary { get; }

		string ToString();
		void BuildSubdomainDictionary();
	}
}
