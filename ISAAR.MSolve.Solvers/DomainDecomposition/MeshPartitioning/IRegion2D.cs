using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning
{
    public enum NodePosition
    {
        Internal, Boundary, External
    }

    public interface IRegion2D
    {
        NodePosition FindRelativePosition(INode node);
    }
}
