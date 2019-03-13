using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public delegate Dictionary<int, SparseVector> NodalLoadsToSubdomainsDistributor(
        Table<INode, DOFType, double> globalNodalLoads);

    public interface IStructuralModel_v2
    {
        Table<INode, DOFType, double> Constraints { get; }
        IReadOnlyList<IElement_v2> Elements { get; }
        IGlobalFreeDofOrdering GlobalDofOrdering { get; set; } //TODO: this should not be managed by the model
        IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; }
        IReadOnlyList<INode> Nodes { get; }
        IReadOnlyList<ISubdomain_v2> Subdomains { get; }

        void AssignLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads); //TODOMaria: Here is where the element loads are assembled
        void AssignMassAccelerationHistoryLoads(int timeStep);
        void ConnectDataStructures();
    }
}
