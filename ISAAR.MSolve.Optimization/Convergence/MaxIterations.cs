using System;

namespace ISAAR.MSolve.Optimization.Convergence
{
    public class MaxIterations : IConvergenceCriterion
    {
        private readonly int maxIterations;

        // Will not be satisfied for iterations 0, 1, .., allowableIterations-1
        public MaxIterations(int maxIterations)
        {
            if (maxIterations < 1)
            {
                throw new ArgumentOutOfRangeException("There must be at least 1 allowable iteration, but there were: "
                                                      + maxIterations);
            }
            this.maxIterations = maxIterations;
        }

        public bool HasConverged(IOptimizationAlgorithm algorithm)
        {
            if (algorithm.CurrentIteration < this.maxIterations)
            {
                return false;
            }
            return true;
        }
    }
}
