using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Bukin's optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class Bukin : SingleObjectiveUnconstrained
    {
        public Bukin()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { -15, -3 };
            this.UpperBound = new double[] { -5, 3 };
            this.ObjectiveFunction = (x) => 100 * Math.Sqrt(Math.Abs(x[1] - 0.01 * Math.Pow(x[0], 2))) +
                    0.01 * Math.Abs(x[0] + 10);
        }
    }
}
