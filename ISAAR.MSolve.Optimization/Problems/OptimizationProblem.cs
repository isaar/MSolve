using System;

namespace ISAAR.MSolve.Optimization.Problems
{
    /// <summary>
    /// Describes the optimization problem to be solved: the objective function(s), the number and bounds of the design 
    /// variables and any constraints. 
    /// </summary>
    public class OptimizationProblem
    {
        /// <summary>
        /// The design factory.
        /// </summary>
        public IDesignFactory DesignFactory
        {
            get; set;
        }

        /// <summary>
        /// The number of continuous (real) design variables.
        /// </summary>
        public int Dimension
        {
            get; set;
        }

        /// <summary>
        /// A vector containing the minimum alloweable values of the continuous (real) design variables. 
        /// Its length must be equal to <see cref="Dimension"/>. 
        /// To represent unbounded design variables, use <see cref="double.MinValue"/>.
        /// </summary>
        public double[] LowerBound
        {
            get; set;
        }

        /// <summary>
        /// A vector containing the maximum alloweable values of the continuous (real) design variables. 
        /// Its length musunbounded design variables, use <see cref="double.MaxValue"/>.
        /// </summary>
        public double[] UpperBound
        {
            get; set;
        }
    }
}