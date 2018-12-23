using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    /// <summary>
    /// This only works if there are no nodal loads or any loading condition other than prescribed displacements.
    /// </summary>
    public class DisplacementControlAnalyzer_v2: NonLinearAnalyzerBase
    {
        private readonly IReadOnlyDictionary<int, IEquivalentLoadsAssembler_v2> equivalentLoadsAssemblers;

        private DisplacementControlAnalyzer_v2(IStructuralModel_v2 model, ISolver_v2 solver, INonLinearProvider_v2 provider,
            IReadOnlyDictionary<int, INonLinearSubdomainUpdater_v2> subdomainUpdaters,
            IReadOnlyDictionary<int, IEquivalentLoadsAssembler_v2> equivalentLoadsAssemblers,
            int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance) :
            base(model, solver, provider, subdomainUpdaters, numIncrements, maxIterationsPerIncrement,
                numIterationsForMatrixRebuild, residualTolerance)
        {
            this.equivalentLoadsAssemblers = equivalentLoadsAssemblers;
        }

        public override void Solve()
        {
            InitializeLogs();

            DateTime start = DateTime.Now;
            UpdateInternalVectors();
            for (int increment = 0; increment < numIncrements; increment++)
            {
                double errorNorm = 0;
                ClearIncrementalSolutionVector();
                UpdateRhs(increment);
                ScaleSubdomainConstraints(increment);

                double firstError = 0;
                int iteration = 0;
                for (iteration = 0; iteration < maxIterationsPerIncrement; iteration++)
                {
                    AddEquivalentNodalLoadsToRHS(increment, iteration);
                    solver.Solve();

                    Dictionary<int, IVector> internalRhsVectors = CalculateInternalRhs(increment, iteration);
                    errorNorm = UpdateResidualForcesAndNorm(increment, internalRhsVectors); // This also sets the rhs vectors in linear systems.
                    if (iteration == 0) firstError = errorNorm;
                    if (errorNorm < residualTolerance) break;
                    //Console.WriteLine($"Increment {increment}, iteration {iteration}: norm2(error) = {errorNorm}");

                    SplitResidualForcesToSubdomains();
                    if ((iteration + 1) % numIterationsForMatrixRebuild == 0)
                    {
                        provider.Reset();
                        BuildMatrices();
                    }
                }
                Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
                SaveMaterialStateAndUpdateSolution();
            }

            // TODO: Logging should be done at each iteration. And it should be done using pull observers
            DateTime end = DateTime.Now;
            StoreLogResults(start, end);
        }

        protected override void InitializeInternalVectors()
        {
            base.InitializeInternalVectors();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                //int idx = FindSubdomainIdx(linearSystems, linearSystem);
                subdomainUpdaters[linearSystem.Subdomain.ID].ScaleConstraints(1 / (double)numIncrements);
            }
        }

        private void AddEquivalentNodalLoadsToRHS(int currentIncrement, int iteration)
        {
            if (iteration != 0)
                return;

            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                //int idx = FindSubdomainIdx(linearSystems, linearSystem);

                double scalingFactor = 1; //((double)currentIncrement + 2) / (currentIncrement + 1); //2; //
                IVector equivalentNodalLoads = equivalentLoadsAssemblers[id].GetEquivalentNodalLoads(u[id], scalingFactor);
                linearSystem.RhsVector.SubtractIntoThis(equivalentNodalLoads);

                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
            }
        }

        // This does nothing at all, as it is written right now
        private void ScaleSubdomainConstraints(int currentIncrement)
        {
            if (currentIncrement == 0)
                return;

            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                //int idx = FindSubdomainIdx(linearSystems, linearSystem);
                double scalingFactor = 1; // ((double)currentIncrement + 2) / (currentIncrement + 1);
                subdomainUpdaters[linearSystem.Subdomain.ID].ScaleConstraints(scalingFactor);
            }
        }

        public class Builder: NonLinearAnalyzerBuilderBase
        {
            private readonly IReadOnlyDictionary<int, IEquivalentLoadsAssembler_v2> equivalentLoadsAssemblers;

            public Builder(IStructuralModel_v2 model, ISolver_v2 solver, INonLinearProvider_v2 provider,
                IReadOnlyDictionary<int, IEquivalentLoadsAssembler_v2> equivalentLoadsAssemblers, int numIncrements):
                base(model, solver, provider, numIncrements)
            {
                MaxIterationsPerIncrement = 1000;
                NumIterationsForMatrixRebuild = 1;
                ResidualTolerance = 1E-3;

                this.equivalentLoadsAssemblers = equivalentLoadsAssemblers;
                //int numSubdomains = model.Subdomains.Count;
                //EquivalentLoadsAssemblers = new EquivalentLoadsAssembler_v2[numSubdomains];
                //for (int i = 0; i < numSubdomains; ++i)
                //{
                //    EquivalentLoadsAssemblers[i] = new EquivalentLoadsAssembler_v2(subdomain, ???); //TODO: ??? must be defined by the provider
                //}
            }

            public DisplacementControlAnalyzer_v2 Build()
            {
                return new DisplacementControlAnalyzer_v2(model, solver, provider, SubdomainUpdaters, equivalentLoadsAssemblers,
                    numIncrements, maxIterationsPerIncrement, numIterationsForMatrixRebuild, residualTolerance);
            }
        }
    }
}
