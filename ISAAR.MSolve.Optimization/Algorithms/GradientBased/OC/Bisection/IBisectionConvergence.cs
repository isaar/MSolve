using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC.Bisection
{
    /// <summary>
    /// Determines if the Optimality Criteria method has converged, based on the value of the lower and upper bisectoning limits.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IBisectionConvergence
    {
        bool HasConverged(double lowerBisectionLimit, double upperBisectionLimit);
    }
}
