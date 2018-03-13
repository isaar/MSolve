using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Holder's table optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class HolderTable : SingleObjectiveUnconstrained
    {
        public HolderTable()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { -10, -10 };
            this.UpperBound = new double[] { 10, 10 };
            this.ObjectiveFunction = (x) => -Math.Abs(Math.Sin(x[0]) * Math.Cos(x[1]) *
                    Math.Exp(Math.Abs(1 - (Math.Sqrt(Math.Pow(x[0], 2) + Math.Pow(x[1], 2)) / Math.PI))));
        }
    }
}