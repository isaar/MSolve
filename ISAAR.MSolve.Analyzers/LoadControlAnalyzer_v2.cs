using System;
using System.Diagnostics;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class LoadControlAnalyzer_v2 : NonLinearAnalyzerBase
    {
        private LoadControlAnalyzer_v2(IStructuralModel_v2 model, ISolver_v2 solver, INonLinearProvider_v2 provider, 
            INonLinearSubdomainUpdater_v2[] subdomainUpdaters,
            int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance) : 
            base(model, solver, provider, subdomainUpdaters, numIncrements, maxIterationsPerIncrement, 
                numIterationsForMatrixRebuild, residualTolerance)
        {
        }

        public override void Solve()
        {
            InitializeLogs();

            DateTime start = DateTime.Now;
            UpdateInternalVectors();//TODOMaria this divides the externally applied load by the number of increments and scatters it to all subdomains and stores it in the class subdomain dictionary and total external load vector
            for (int increment = 0; increment < numIncrements; increment++)
            {
                double errorNorm = 0;
                ClearIncrementalSolutionVector();//TODOMaria this sets du to 0
                UpdateRhs(increment);//TODOMaria this copies the residuals stored in the class dictionary to the subdomains

                double firstError = 0;
                int iteration = 0;
                for (iteration = 0; iteration < maxIterationsPerIncrement; iteration++)
                {
                    solver.Solve();
                    errorNorm = globalRhsNorm != 0 ? CalculateInternalRhs(increment, iteration) / globalRhsNorm : 0;// (rhsNorm*increment/increments) : 0;//TODOMaria this calculates the internal force vector and subtracts it from the external one (calculates the residual)
                    Console.WriteLine($"Increment {increment}, iteration {iteration}: norm2(error) = {errorNorm}");

                    if (iteration == 0) firstError = errorNorm;

                    if (IncrementalDisplacementsLog != null) IncrementalDisplacementsLog.StoreDisplacements_v2(uPlusdu);

                    if (errorNorm < residualTolerance) break;

                    SplitResidualForcesToSubdomains();//TODOMaria scatter residuals to subdomains
                    if ((iteration + 1) % numIterationsForMatrixRebuild == 0) // Matrix rebuilding should be handled in another way. E.g. in modified NR, it must be done at each increment.
                    {
                        provider.Reset();
                        BuildMatrices();
                    }
                }
                Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
                SaveMaterialStateAndUpdateSolution();
            }
            //            ClearMaterialStresses();

            // TODO: Logging should be done at each iteration. And it should be done using pull observers
            DateTime end = DateTime.Now;
            StoreLogResults(start, end);
        }

        public class Builder: NonLinearAnalyzerBuilderBase
        {
            public Builder(IStructuralModel_v2 model, ISolver_v2 solver, INonLinearProvider_v2 provider, int numIncrements):
                base(model, solver, provider, numIncrements)
            {
                MaxIterationsPerIncrement = 1000;
                NumIterationsForMatrixRebuild = 1;
                ResidualTolerance = 1E-8;
            }

            public LoadControlAnalyzer_v2 Build()
            {
                return new LoadControlAnalyzer_v2(model, solver, provider, SubdomainUpdaters,
                    numIncrements, maxIterationsPerIncrement, numIterationsForMatrixRebuild, residualTolerance);
            }
        }
    }
}
