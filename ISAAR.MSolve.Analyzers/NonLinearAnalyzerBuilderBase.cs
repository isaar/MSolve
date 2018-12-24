using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public abstract class NonLinearAnalyzerBuilderBase
    {
        protected int maxIterationsPerIncrement = 1000;
        protected readonly IStructuralModel_v2 model;
        protected readonly int numIncrements;
        protected int numIterationsForMatrixRebuild = 1;
        protected readonly INonLinearProvider_v2 provider;
        protected double residualTolerance = 1e-8;
        protected readonly ISolver_v2 solver;

        protected NonLinearAnalyzerBuilderBase(IStructuralModel_v2 model, ISolver_v2 solver, INonLinearProvider_v2 provider, 
            int numIncrements)
        {
            //TODO: this should belong to all (child) analyzer builders
            this.model = model;
            this.solver = solver;
            this.provider = provider;
            this.numIncrements = numIncrements;
            SubdomainUpdaters = CreateDefaultSubdomainUpdaters();
        }

        public int MaxIterationsPerIncrement
        {
            get => maxIterationsPerIncrement;
            set
            {
                if (value < 1) throw new ArgumentException("Max iterations per increment must be >= 1");
                maxIterationsPerIncrement = value;
            }
        }

        public int NumIterationsForMatrixRebuild
        {
            get => numIterationsForMatrixRebuild;
            set
            {
                if (value < 1) throw new ArgumentException("Iterations number for matrix rebuild must be >= 1");
                numIterationsForMatrixRebuild = value;
            }
        }

        public double ResidualTolerance
        {
            get => residualTolerance;
            set
            {
                if (value <= 0.0) throw new ArgumentException("Residual tolerance must be positive");
                residualTolerance = value;
            }
        }

        public IReadOnlyDictionary<int, INonLinearSubdomainUpdater_v2> SubdomainUpdaters { get; set; }

        private IReadOnlyDictionary<int, INonLinearSubdomainUpdater_v2> CreateDefaultSubdomainUpdaters()
        {
            int numSubdomains = model.Subdomains.Count;
            var subdomainUpdaters = new Dictionary<int, INonLinearSubdomainUpdater_v2>(numSubdomains);
            for (int i = 0; i < numSubdomains; ++i)
            {
                subdomainUpdaters[i] = new NonLinearSubdomainUpdater_v2(model.Subdomains[i]);
            }
            return subdomainUpdaters;

        }
    }
}
