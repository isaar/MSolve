using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Levi's No.13 optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class Levy13 : SingleObjectiveUnconstrained
    {
        public Levy13()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { -10, -10 };
            this.UpperBound = new double[] { 10, 10 };
            this.ObjectiveFunction = (x) => Math.Pow(Math.Sin(3 * Math.PI * x[0]), 2) +
                    Math.Pow(x[0] - 1, 2) * (1 + Math.Pow(Math.Sin(3 * Math.PI * x[1]), 2)) +
                    Math.Pow(x[1] - 1, 2) * (1 + Math.Pow(Math.Sin(2 * Math.PI * x[1]), 2));
        }
    }
}
