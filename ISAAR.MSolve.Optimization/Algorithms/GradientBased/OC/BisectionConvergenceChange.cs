using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC
{
    /// <summary>
    /// OC convergences when (upper-lower) &lt;= tol. Presented in "A 99 line topology optimization code written in Matlab - 
    /// O. Sigmund - 1991".
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class BisectionConvergenceChange : IBisectionConvergence
    {
        private readonly double tolerance;

        public BisectionConvergenceChange(double tolerance) => this.tolerance = tolerance;

        public bool HasConverged(double lowerBisectionLimit, double upperBisectionLimit)
            => upperBisectionLimit - lowerBisectionLimit <= tolerance;
    }
}
