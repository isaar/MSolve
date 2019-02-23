using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Analyzers.Multiscale
{
    public interface IKinematicRelationsStrategy
    {
        IMatrixView GetKinematicRelationsMatrix();
    }
}
