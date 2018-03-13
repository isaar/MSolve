using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the McCormick's optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class McCormick : SingleObjectiveUnconstrained
    {
        public McCormick()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { -1.5, -3.0 };
            this.UpperBound = new double[] { 4.0, 4.0 };
            this.ObjectiveFunction = (x) => Math.Sin(x[0] + x[1]) + Math.Pow(x[0] - x[1], 2) - 1.5 * x[0] + 2.5 * x[1] + 1;
        }
    }
}