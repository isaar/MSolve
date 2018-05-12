using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IStructuralModel
    {
		Dictionary<int, ISubdomain> ISubdomainsDictionary { get;  }
	    IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; }
		void AssignLoads();
		void AssignMassAccelerationHistoryLoads(int timeStep);
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
