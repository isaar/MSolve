using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.Commons;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IStructuralModel_v2
    {
        Table<INode, DOFType, double> Constraints { get; }
        IReadOnlyList<IElement> Elements { get; }
        IGlobalFreeDofOrdering GlobalDofOrdering { get; } //TODO: this should not be managed by the model
        IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; }
        IReadOnlyList<INode> Nodes { get; }
        //Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary { get; } 
        IReadOnlyList<ISubdomain_v2> Subdomains { get; }

        void AssignLoads();
        void AssignMassAccelerationHistoryLoads(int timeStep);
	}
}
