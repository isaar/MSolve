using System;
using ISAAR.MSolve.Optimization.Problems;
using System.Linq;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Styblinski-Tang's optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class StyblinskiTang : SingleObjectiveUnconstrained
    {
        public StyblinskiTang()
        {
            this.Dimension = 10;

            this.LowerBound = new double[Dimension];
            LowerBound = LowerBound.Select(i => -10.0).ToArray();

            this.UpperBound = new double[Dimension];
            UpperBound = UpperBound.Select(i => 10.0).ToArray();

            this.ObjectiveFunction = (x) =>
            {
                double f = 0.0;
                for (int i = 0; i < x.Length; i++)
                {
                    f += Math.Pow(x[i], 4) - 16 * Math.Pow(x[i], 2) + 5 * x[i];
                }
                return f / 2.0;
            };
        }
    }
}