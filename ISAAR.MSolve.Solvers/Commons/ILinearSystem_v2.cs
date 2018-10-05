using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

//TODO: perhaps this, the providers and analyzers should be generic.
//TODO: This should hide distributed systems arising in domain decomposition methods. That should be transparent to the analyzers
//      and providers.
namespace ISAAR.MSolve.Solvers.Commons
{
    public interface ILinearSystem_v2
    {
        int ID { get; } //TODO: delete this once subdomains have been abstracted.

        //TODO: this is error prone. The implementation should manage the state, by funneling access to the matrix.
        bool IsMatrixModified { get; set; }

        IMatrix Matrix { get; }
        IVector RhsVector { get; }
        IVector Solution { get; }
    }
}
