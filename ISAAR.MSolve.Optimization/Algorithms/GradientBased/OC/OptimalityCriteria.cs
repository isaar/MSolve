using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Problems;

//TODO: Use and enrich the auxilliary classes for other optimization algorithms.
namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC
{
    /// <summary>
    /// Implements the Optimality Criteria optimization algorithm for problems with a differentiable objective function and a
    /// single monotonic equality constraint function.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class OptimalityCriteria
    {
        private readonly IBisectionConvergence bisectionConvergence;
        private readonly double bisectionInitialLimitLower, bisectionInitialLimitUpper;
        private readonly double boundLower, boundUpper;
        private readonly EqualityConstraint constraint;
        private readonly double dampingCoeff;
        private readonly double moveLimit;
        private readonly DifferentiableObjectiveFunction objective;

        internal OptimalityCriteria(DifferentiableObjectiveFunction objective, EqualityConstraint constraint,
            double boundLower, double boundUpper, double bisectionInitialLimitLower, double bisectionInitialLimitUpper,
            IBisectionConvergence bisectionConvergence, double dampingCoeff, double moveLimit)
        {
            this.objective = objective;
            this.constraint = constraint;
            this.boundLower = boundLower;
            this.boundUpper = boundUpper;
            this.bisectionInitialLimitLower = bisectionInitialLimitLower;
            this.bisectionInitialLimitUpper = bisectionInitialLimitUpper;
            this.dampingCoeff = dampingCoeff;
            this.moveLimit = moveLimit;
            this.bisectionConvergence = bisectionConvergence;
        }

        public (double fMin, Vector xBest) Optimize(Vector xInitial)
        {
            Vector x = xInitial.Copy();
            Vector xNext;
            double f = double.NaN;
            Vector gradF = null;
            double xChange = 1.0;
            while (xChange > 0.01)
            {
                (f, gradF) = objective(x);
                xNext = UpdateDesign(x, gradF);
                xChange = (xNext - x).MaxAbsolute();
                x = xNext;
            }
            return (f, x);
        }

        //TODO: what about passive elements? They should not be design variables at all.
        private Vector UpdateDesign(Vector x, Vector gradF) 
        {
            double l1 = bisectionInitialLimitLower;
            double l2 = bisectionInitialLimitUpper;
            Vector xNext = null;
            while (l2 - l1 > 1E-4)
            {
                double lmid = 0.5 * (l1 + l2);

                // Update densities
                xNext = x.DoEntrywise(gradF, (xi, gradFi) =>
                {
                    double xiBi = xi * Math.Pow(-gradFi / lmid, dampingCoeff);
                    double xiLow = Math.Max(boundLower, xi - moveLimit);
                    double xiHigh = Math.Min(boundUpper, xi + moveLimit);
                    if (xiBi <= xiLow) return xiLow;
                    else if (xiBi < xiHigh) return xiBi;
                    else return xiHigh;
                });

                // Bi-sectioning
                if (constraint(xNext) > 0) l1 = lmid;
                else l2 = lmid;
            }
            return xNext;
        }
    }
}
