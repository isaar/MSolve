using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling
{
    public class ExponentialFitnessScaling<T> : IFitnessScalingStrategy<T>
    {
        // The higher it is the more biased towards the fittest. 0 means equal expectations. Negative means bias towards least fit.
        private readonly double biasTowardFittest;

        public ExponentialFitnessScaling(double biasTowardFittest)
        {
            this.biasTowardFittest = biasTowardFittest;
        }

        public double[] CalculateExpectations(Individual<T>[] individuals)
        {
            // Find max fitness
            double maxFitness = individuals[0].Fitness;
            foreach (var individual in individuals) 
            {
                if (individual.Fitness < maxFitness) maxFitness = individual.Fitness;
            }

            // Assign expectations
            double[] expectations = new double[individuals.Length];
            for (int i = 0; i < individuals.Length; ++i)
            {
                expectations[i] = Math.Exp(-biasTowardFittest * individuals[i].Fitness / maxFitness);
            }
            return expectations;
        }
    }
}
