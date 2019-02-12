using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Problems;

//TODO: Use and enrich the auxilliary classes for other optimization algorithms.
namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC
{
    //TODO: use a design class for these.
    public delegate (double f, Vector gradF) DifferentiableObjectiveFunction(Vector x);
    public delegate double EqualityConstraint(Vector x);

    /// <summary>
    /// Implements the Optimality Criteria optimization algorithm for problems with a differentiable objective function and a
    /// single monotonic equality constraint function.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class OptimalityCriteria
    {
        private readonly DifferentiableObjectiveFunction objective;
        private readonly EqualityConstraint constraint;
        private readonly double lowerBound;
        private readonly double upperBound;

        public OptimalityCriteria(DifferentiableObjectiveFunction objective, EqualityConstraint constraint,
            double lowerBound, double upperBound)
        {
            this.objective = objective;
            this.constraint = constraint;
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
        }

        public double DampingCoeff { get; set; } = 0.5;
        public double MaxDesignVariableChange { get; set; } = 0.01;
        public double MoveLimit { get; set; } = 0.2;

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
            double l1 = 0.0;
            double l2 = 1E5;
            Vector xNext = null;
            while (l2 - l1 > 1E-4)
            {
                double lmid = 0.5 * (l1 + l2);

                // Update densities
                xNext = x.DoEntrywise(gradF, (xi, gradFi) =>
                {
                    double xiBi = xi * Math.Pow(-gradFi / lmid, DampingCoeff);
                    double xiLow = Math.Max(lowerBound, xi - MoveLimit);
                    double xiHigh = Math.Min(upperBound, xi + MoveLimit);
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
