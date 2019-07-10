using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC.Convergence
{
    public class MaxIterationsConvergence : IOptimalityCriteriaConvergence
    {
        private readonly int maxIterations;

        /// <summary>
        /// </summary>
        /// <param name="maxIterations">
        /// This convergence criterion is not satisfied if 0 &lt;= currentIteration &lt; <paramref name="maxIterations"/>.
        /// </param>
        public MaxIterationsConvergence(int maxIterations)
        {
            this.maxIterations = maxIterations;
        }

        public bool HasConverged(int currentIteration, double currentObjectiveFunction, IVectorView nextDesignVariables)
            => currentIteration >= maxIterations; // Other criteria may cause the > case.
    }
}
