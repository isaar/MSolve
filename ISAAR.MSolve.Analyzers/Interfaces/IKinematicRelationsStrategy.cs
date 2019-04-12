using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IKinematicRelationsStrategy
    {
        double[,] GetNodalKinematicRelationsMatrix(INode boundaryNode);
    }
}
