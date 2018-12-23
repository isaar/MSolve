using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISolver_v2 : ISystemMatrixObserver
    {
        IDofOrderer DofOrderer { get; }

        IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        //TODO: Ideally the provider/analyzer will not even have to pass the subdomain.
        IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider elementMatrixProvider);

        // This needs to be called only once (by the analyzer), not every time the matrix must be factorized/modified, since it  
        // could potentially perform actions that must not be repeated or are too expesive. Instead modifying the system matrix, 
        // should also notify the solver.
        void Initialize();
        void Solve();
    }
}
