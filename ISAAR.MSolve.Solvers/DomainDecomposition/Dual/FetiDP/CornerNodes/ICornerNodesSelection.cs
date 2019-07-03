using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes
{
    public interface ICornerNodeSelection
    {
        //TODO: These should probably be HashSet instead of array
        Dictionary<int, HashSet<INode>> SelectCornerNodesOfSubdomains();
    }
}
