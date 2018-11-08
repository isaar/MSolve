//using System.Collections.Generic;
//using ISAAR.MSolve.Discretization.Interfaces;
//using ISAAR.MSolve.LinearAlgebra.Matrices;
//using ISAAR.MSolve.Solvers.Commons;

//namespace ISAAR.MSolve.Solvers.Interfaces
//{
//    public interface ISolver_v2
//    {
//        IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

//        //TODO: Ideally the provider/analyzer will not even have to pass the subdomain.
//        IMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider elementMatrixProvider);

//        //TODO: I think this needs to be called only once (by the analyzer), not every time the matrix must be factorized.
//        void Initialize();

//        void Solve();
//    }
//}
