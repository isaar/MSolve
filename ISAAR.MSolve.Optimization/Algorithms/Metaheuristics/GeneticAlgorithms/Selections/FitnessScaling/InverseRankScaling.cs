using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling
{
    public class InverseRankScaling<T> : IFitnessScalingStrategy<T>
    {
        // The higher it is the more biased towards the fittest. 0 means equal expectations. Negative means bias towards least fit.
        private readonly double exponent;

        public InverseRankScaling(double exponent = 0.5)
        {
            this.exponent = exponent;
        }

        public double[] CalculateExpectations(Individual<T>[] individuals)
        {
            Array.Sort(individuals);
            double[] expectations = new double[individuals.Length];
            for (int i = 0; i < individuals.Length; ++i)
            {
                expectations[i] = 1.0 / Math.Pow(i + 1, exponent);
            }
            return expectations;
        }
    }
}
