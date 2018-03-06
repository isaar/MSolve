using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Cross-In-Tray's optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class CrossInTray : SingleObjectiveUnconstrained
    {
        public CrossInTray()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { -10, -10 };
            this.UpperBound = new double[] { 10, 10 };
            this.ObjectiveFunction = (x) => -0.0001 * Math.Pow(Math.Abs(Math.Sin(x[0]) * Math.Sin(x[1]) *
                    Math.Exp(Math.Abs(100 - (Math.Sqrt(Math.Pow(x[0], 2) + Math.Pow(x[1], 2))) / Math.PI))) + 1, 0.1);
        }
    }
}
