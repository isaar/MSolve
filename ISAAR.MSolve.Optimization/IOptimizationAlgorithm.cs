using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization
{
    public interface IOptimizationAlgorithm
    {
        double BestFitness
        {
            get;
        }

        double[] BestPosition
        {
            get;
        }

        int CurrentIteration
        { 
            get; 
        }

        double CurrentFunctionEvaluations
        {
            get;
        }

        void Solve();
    }
}