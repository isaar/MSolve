using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Matyas's optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class Matyas : SingleObjectiveUnconstrained
    {
        public Matyas()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { -10, -10 };
            this.UpperBound = new double[] { 10, 10 };
            this.ObjectiveFunction = (x) => 0.26 * (Math.Pow(x[0], 2) + Math.Pow(x[1], 2)) + 0.48 * x[0] * x[1];
        }
    }
}
