using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: merge this code with the IConvergenceCriterion defined for general algorithms.
//TODO: not sure if the designVariables vector must correspond to the next iteration, while the fitness and iteration are current
namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC.Convergence
{
    public interface IOptimalityCriteriaConvergence
    {
        /// <summary>
        /// </summary>
        /// <param name="iteration">
        /// The current iteration, starting from 0. At each iteration the objective function is evaluated. Usually iteration 0
        /// is initialization, although this may change if the objective function is not evaluated.
        /// </param>
        /// <param name="designVariables"></param>
        /// <param name="objectiveFunction"></param>
        bool HasConverged(int currentIteration, double currentObjectiveFunction, IVectorView nextDesignVariables);
    }
}
