using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC.Bisection
{
    /// <summary>
    /// OC convergences when (upper-lower)/(upper+lower) &lt;= tol. Presented in "A 99 line topology optimization code written   
    /// in Matlab - O. Sigmund - 1991".
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class BisectionConvergenceRelativeChange : IBisectionConvergence
    {
        private readonly double tolerance;

        public BisectionConvergenceRelativeChange(double tolerance) => this.tolerance = tolerance;

        public bool HasConverged(double lowerBisectionLimit, double upperBisectionLimit)
            => (upperBisectionLimit - lowerBisectionLimit) / (upperBisectionLimit + lowerBisectionLimit) <= tolerance;
    }
}
