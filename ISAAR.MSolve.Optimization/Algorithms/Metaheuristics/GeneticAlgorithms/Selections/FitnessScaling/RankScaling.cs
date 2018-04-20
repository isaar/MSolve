using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling
{
    public class RankScaling<T> : IFitnessScalingStrategy<T>
    {
        // The higher it is the more biased towards the fittest. 0 means equal expectations. Negative means bias towards least fit.
        private readonly double exponent;

        public RankScaling(double exponent = 1.0)
        {
            this.exponent = exponent;
        }

        public double[] CalculateExpectations(Individual<T>[] individuals)
        {
            Array.Sort(individuals);
            int count = individuals.Length;
            double[] expectations = new double[count];
            for (int i = 0; i < count; ++i)
            {
                expectations[i] = Math.Pow(count - i, exponent); // Fitest rank = N^p, Least fit rank = [N-(N-1)]^p = 1
            }
            return expectations;
        }
    }
}
