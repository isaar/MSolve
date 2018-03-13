using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling
{
    public interface IFitnessScalingStrategy<T>
    {
        /// <summary>
        /// Computes the expectation of each <see cref="Individual{T}"/>. They must be positive numbers. 
        /// These expectations will serve as the probability that each <see cref="Individual{T}"/> 
        /// will be selected for reproduction.
        /// </summary>
        /// <param name="individuals"></param>
        /// <returns>A vector with the expectations</returns>
        double[] CalculateExpectations(Individual<T>[] individuals);
    }
}
