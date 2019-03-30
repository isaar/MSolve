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
    public class StaticAnalyzer_v2 : INonLinearParentAnalyzer_v2
    {
        private readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        private readonly IStructuralModel_v2 model;
        private readonly IStaticProvider_v2 provider;
        private readonly ISolver_v2 solver;

        public StaticAnalyzer_v2(IStructuralModel_v2 model, ISolver_v2 solver, IStaticProvider_v2 provider, 
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
            if (isFirstAnalysis)
            {
                // The order in which the next initializations happen is very important.
                model.ConnectDataStructures();
                solver.OrderDofs(false);
                foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
                {
                    linearSystem.Reset(); // Necessary to define the linear system's size 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
            }
            else
            {
                foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
                {
                    //TODO: Perhaps these shouldn't be done if an analysis has already been executed. The model will not be 
                    //      modified. Why should the linear system be?
                    linearSystem.Reset(); 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
            }

            //TODO: Perhaps this should be called by the child analyzer
            BuildMatrices(); 

            // Loads must be created after building the matrices.
            //TODO: Some loads may not have to be recalculated each time the stiffness changes.
            model.AssignLoads(solver.DistributeNodalLoads); 
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }

            if (ChildAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            ChildAnalyzer.Initialize(isFirstAnalysis);
        }

        public void Solve()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            ChildAnalyzer.Solve();
        }
    }
}
