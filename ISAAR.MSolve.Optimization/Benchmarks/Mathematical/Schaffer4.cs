using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Schaffer's No.4 optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class Schaffer4 : SingleObjectiveUnconstrained
    {
        public Schaffer4()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { -100, -100 };
            this.UpperBound = new double[] { 100, 100 };
            this.ObjectiveFunction = (x) => 0.5 
                + (Math.Pow(Math.Cos(Math.Sin(Math.Abs(Math.Pow(x[0], 2) - Math.Pow(x[1], 2)))), 2) - 0.5)
                / Math.Pow(1 + 0.001 * (Math.Pow(x[0], 2) + Math.Pow(x[1], 2)), 2);
        }
    }
}
