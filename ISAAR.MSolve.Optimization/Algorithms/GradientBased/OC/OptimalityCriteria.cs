using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC.Bisection;
using ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC.Convergence;
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
        private readonly IOptimalityCriteriaConvergence ocConvergence;
        private readonly DifferentiableObjectiveFunction objective;

        internal OptimalityCriteria(DifferentiableObjectiveFunction objective, EqualityConstraint constraint,
            double boundLower, double boundUpper, IOptimalityCriteriaConvergence ocConvergence, 
            double bisectionInitialLimitLower, double bisectionInitialLimitUpper, IBisectionConvergence bisectionConvergence,
            double dampingCoeff, double moveLimit)
        {
            this.objective = objective;
            this.constraint = constraint;
            this.boundLower = boundLower;
            this.boundUpper = boundUpper;
            this.ocConvergence = ocConvergence;
            this.bisectionInitialLimitLower = bisectionInitialLimitLower;
            this.bisectionInitialLimitUpper = bisectionInitialLimitUpper;
            this.dampingCoeff = dampingCoeff;
            this.moveLimit = moveLimit;
            this.bisectionConvergence = bisectionConvergence;
        }

        public (Vector xBest, double fMin) Optimize(Vector xInitial)
        {
            Vector x = xInitial.Copy(); //TODO: is this necessary?
            double f = double.NaN;
            Vector gradF = null;
            int iteration = -1;
            do
            {
                ++iteration;
                (f, gradF) = objective(x);
                x = UpdateDesign(x, gradF);
            } while (!ocConvergence.HasConverged(iteration, f, x));
            return (x, f);
        }

        //TODO: what about passive elements? They should not be design variables at all.
        private Vector UpdateDesign(Vector x, Vector gradF) 
        {
            double l1 = bisectionInitialLimitLower;
            double l2 = bisectionInitialLimitUpper;
            Vector xNext = null;
            while (!bisectionConvergence.HasConverged(l1, l2))
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
