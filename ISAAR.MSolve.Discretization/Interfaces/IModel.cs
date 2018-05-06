using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IModel
    {
		Dictionary<int,ISubdomain> ISubdomainsDictionary { get; }
		Dictionary<int, INode> NodesDictionary { get; }
		Dictionary<int, IElement> ElementsDictionary { get;}
		Dictionary<int, ISubdomain> SubdomainsDictionary { get;  }
		//Dictionary<int, Cluster> ClustersDictionary { get; set; }

		IList<INode> Nodes { get; }
		IList<IElement> Elements { get; }
		IList<ISubdomain> Subdomains { get; }
		Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary { get; }
		int TotalDOFs { get; }

		void AssignLoads();
		void AssignMassAccelerationHistoryLoads();
		void ConnectDataStructures();
		void Clear();

	}

	public enum DOFType
	{
		Unknown = 0,
		X = 1,
		Y = 2,
		Z = 3,
		RotX = 4,
		RotY = 5,
		RotZ = 6,
		Pore = 7
	}

}
