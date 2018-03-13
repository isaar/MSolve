using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Problems
{
    /// <summary>
    /// Utility class that handles all logic concerning the validity of the <see cref="OptimizationProblem"/>'s properties 
    /// inputted by the user. 
    /// </summary>
    static class ProblemChecker
    {
        /// <summary>
        /// Checks all properties of the provided <see cref="OptimizationProblem"/>. If any violation is found, an
        /// <see cref="ArgumentException"/> will be thrown.
        /// </summary>
        /// <param name="problem">The <see cref="OptimizationProblem"/> whose properties will be checked.</param>
        /// <exception cref="ArgumentException">Thrown when a property's value is not legal 
        ///     or the values of 2 or more properties are legal, but cannot be used together.</exception>
        public static void Check(OptimizationProblem problem)
        {
            CheckDesignVariables(problem);
            CheckBounds(problem);
            CheckObjectives(problem);
        }

        /// <summary>
        /// Checks the number of design variables and the length of the bound constraint vectors.
        /// </summary>
        /// <param name="problem">The <see cref="OptimizationProblem"/> whose properties will be checked.</param>
        private static void CheckDesignVariables(OptimizationProblem problem)
        {
            if (problem.Dimension < 1)
            {
                throw new ArgumentException("The number of continuous design variables must be >= 1, but was : " 
                    + problem.Dimension);
            }
            if (problem.LowerBound.Length != problem.Dimension)
            {
                throw new ArgumentException("There number of continuous lower bounds was " + problem.LowerBound.Length + 
                    " , which is different than the number of continuous design variables " + problem.Dimension 
                    + ". They must be the same");
            }
            if (problem.UpperBound.Length != problem.Dimension)
            {
                throw new ArgumentException("There number of continuous upper bounds was " + problem.LowerBound.Length +
                    " , which is different than the number of continuous design variables " + problem.Dimension
                    + ". They must be the same");
            }
        }

        /// <summary>
        /// Checks the contents of the bound constraint vectors.
        /// </summary>
        /// <param name="problem">The <see cref="OptimizationProblem"/> whose properties will be checked.</param>
        private static void CheckBounds(OptimizationProblem problem)
        {
            for (int i = 0; i < problem.Dimension; i++)
            {
                if (problem.LowerBound[i] >= problem.UpperBound[i])
                {
                    throw new ArgumentException("Lower bound value of design variable " + i
                        + " must be lower than the corresponding upper bound!");
                }
            }
        }

        /// <summary>
        /// Checks the objective functions.
        /// </summary>
        /// <param name="problem">The <see cref="OptimizationProblem"/> whose properties will be checked.</param>
        private static void CheckObjectives(OptimizationProblem problem)
        {
           
        }
    }
}
