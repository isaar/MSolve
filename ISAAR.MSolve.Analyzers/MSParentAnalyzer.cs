using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Analyzers
{
    /// <summary>
    /// Parent Analyzer for nonlinear static problems assosiated with a microstructure bvp handled by MicrostructureBvpNRNLAnalyzer
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class MSParentAnalyzer : INonLinearParentAnalyzer_v2
    {
        private readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        private readonly IStructuralModel_v2 model;
        private readonly IStaticProvider_v2 provider;
        private readonly ISolver_v2 solver;
        private bool firstInitialization;

        public MSParentAnalyzer(IStructuralModel_v2 model, ISolver_v2 solver, IStaticProvider_v2 provider, 
            IChildAnalyzer childAnalyzer)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.provider = provider;
            this.ChildAnalyzer = childAnalyzer;
            this.ChildAnalyzer.ParentAnalyzer = this;
        }

        public Dictionary<int, IAnalyzerLog_v2[]> Logs { get; } = new Dictionary<int, IAnalyzerLog_v2[]>();

        public IChildAnalyzer ChildAnalyzer { get; }

        public void BuildMatrices()
        {
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                linearSystem.Matrix = provider.CalculateMatrix(linearSystem.Subdomain);
            }
        }

        public IVector GetOtherRhsComponents(ILinearSystem_v2 linearSystem, IVector currentSolution)
        {
            //TODO: use a ZeroVector class that avoid doing useless operations or refactor this method. E.g. let this method 
            // alter the child analyzer's rhs vector, instead of the opposite (which is currently done).
            return linearSystem.CreateZeroVector();
        }

        public void Initialize(bool isFirstAnalysis = true)
        {

            if (ChildAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            //InitalizeMatrices();

            //model.ConnectDataStructures();
            //solver.OrderDofsAndClearLinearSystems(); //TODO MS Since orderdofs is called here it cannot be called twice (hence in subdomainmicrobase) so modify
            //solver.ResetSubdomainForcesVector();

            //model.AssignLoads(); //mhdenizei ta subdomain.forces ean exoume epivalei fortio afou tous pernaei ta fortia

            //TODO: this should be done elsewhere. It makes sense to assign the RHS vector when the stiffness matrix is assigned
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }

            ChildAnalyzer.Initialize(isFirstAnalysis);
        }

        public void Solve()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            BuildMatrices(); //TODO: this should be called by the class that calls model.AssignLoads() and before it. 
            ChildAnalyzer.Solve();
        }
    }
}
