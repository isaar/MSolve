using System;
using ISAAR.MSolve.Optimization.Problems;
using System.Linq;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Rosenbrock's optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class Rosenbrock : SingleObjectiveUnconstrained
    {
        public Rosenbrock()
        {
            this.Dimension = 10;

            this.LowerBound = new double[Dimension];
            LowerBound = LowerBound.Select(i => -10.0).ToArray();

            this.UpperBound = new double[Dimension];
            UpperBound = UpperBound.Select(i => 10.0).ToArray();

            this.ObjectiveFunction = (x) =>
            {
                double f = 0.0;
                double t1;
                double t2;

                for (int i = 0; i < (x.Length - 1); i++)
                {
                    t1 = (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]);
                    t2 = (1 - x[i]) * (1 - x[i]);
                    f += 100.0 * t1 + t2;
                }

                return f;
            };
        }
    }
}