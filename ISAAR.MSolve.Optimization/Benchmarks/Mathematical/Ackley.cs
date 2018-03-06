using System;
using System.Linq;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Ackley's optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class Ackley : SingleObjectiveUnconstrained
    {
        public Ackley(int dim)
        {
            this.Dimension = dim;

            this.LowerBound = new double[Dimension];
            LowerBound = LowerBound.Select(i => -5.0).ToArray();

            this.UpperBound = new double[Dimension];
            UpperBound = UpperBound.Select(i => 5.0).ToArray();

            this.ObjectiveFunction = (x) =>
            {
                double rootSum = 0.0;
                for (int i = 0; i < this.Dimension; ++i)
                {
                    rootSum += Math.Pow(x[i], 2);
                }
                rootSum /= this.Dimension;

                double expSum = 0.0;
                for (int i = 0; i < this.Dimension; ++i)
                {
                    expSum += Math.Cos(2 * Math.PI * x[i]);
                }
                expSum /= this.Dimension;

                return -20 * Math.Exp(-0.2 * Math.Sqrt(rootSum)) - Math.Exp(expSum) + Math.E + 20;
            };
        }
    }
}