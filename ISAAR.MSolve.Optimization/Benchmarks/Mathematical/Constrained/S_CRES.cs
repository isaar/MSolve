using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical.Constrained
{
    public class S_CRES : SingleObjectiveConstrained
    {
        public S_CRES()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { 0, 0 };
            this.UpperBound = new double[] { 6, 6 };

            ObjectiveFunction = (x) =>
            {
               double x1 = x[0];
               double x2 = x[1];
               double x1_2 = Math.Pow(x1, 2.0);
               double x2_2 = Math.Pow(x2, 2.0);
               return Math.Pow(x1_2 + x2 - 11, 2) + Math.Pow(x1 + x2_2 - 7, 2);
            };

            AddConstraintFunction((x) => Math.Pow(x[0] - 0.05, 2) + Math.Pow(x[1] - 2.5, 2) - 4.84);
            AddConstraintFunction((x) => -Math.Pow(x[0], 2) - Math.Pow(x[1] - 2.5, 2) + 4.84);
        }
    }
}