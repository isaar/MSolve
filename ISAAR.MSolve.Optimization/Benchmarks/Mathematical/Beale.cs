using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Beale's optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class Beale : SingleObjectiveUnconstrained
    {
        public Beale()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { -4.5, -4.5 };
            this.UpperBound = new double[] { 4.5, 4.5 };
            this.ObjectiveFunction = (x) => Math.Pow((1.5 - x[0] + x[0] * x[1]), 2)
                    + Math.Pow((2.25 - x[0] + x[0] * x[1] * x[1]), 2)
                    + Math.Pow((2.625 - x[0] + x[0] * x[1] * x[1] * x[1]), 2);
        }
    }
}
