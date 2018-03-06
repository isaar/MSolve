using System;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Mathematical
{
    /// <summary>
    /// Class for the Booth's optimization problem.
    /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization">Wikipedia: Test functions for optimization</see>
    /// </summary>
    public class Booth : SingleObjectiveUnconstrained
    {
        public Booth()
        {
            this.Dimension = 2;
            this.LowerBound = new double[] { -10, -10 };
            this.UpperBound = new double[] { 10, 10 };
            this.ObjectiveFunction = (x) => 
                Math.Pow(x[0] + 2 * x[1] - 7, 2) + Math.Pow(2 * x[0] + x[1] - 5, 2);
        }
    }
}