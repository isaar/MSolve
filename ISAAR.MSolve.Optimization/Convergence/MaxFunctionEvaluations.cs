using System;

namespace ISAAR.MSolve.Optimization.Convergence
{
    public class MaxFunctionEvaluations : IConvergenceCriterion
    {
        private readonly int maxFES;

        public MaxFunctionEvaluations(int maxFES)
        {
            //if (maxFES < 100000*dimension)
            //{
            //    throw new ArgumentOutOfRangeException("There must be at least one function evaluation, but there were: "
            //                                          + maxFES);
            //}
            this.maxFES = maxFES;
        }

        public bool HasConverged(IOptimizationAlgorithm algorithm)
        {
            if (algorithm.CurrentFunctionEvaluations < this.maxFES)
            {
                return false;
            }
            return true;
        }
    }
}