using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling;
using ISAAR.MSolve.Optimization.Commons;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections
{
    public abstract class AbstractRemainderSelection<T> : ISelectionStrategy<T>
    {
        private readonly IFitnessScalingStrategy<T> fitnessScaling;
        protected IGenerator RNG { get; }

        protected AbstractRemainderSelection(IFitnessScalingStrategy<T> fitnessScaling,
            IGenerator randomNumberGenerator)
        {
            if (fitnessScaling == null) throw new ArgumentException("The fitness scaling strategy must not be null");
            this.fitnessScaling = fitnessScaling;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.RNG = randomNumberGenerator;
        }

        public Individual<T>[][] Apply(Individual<T>[] population, int parentGroupsCount, int parentsPerGroup)
        {
            int totalParentsCount = parentGroupsCount * parentsPerGroup;
            var allParents = new Individual<T>[totalParentsCount];

            // Integral phase
            double[] fractionalParts;
            LinkedList<Individual<T>> someParents = IntegralPhase(population, totalParentsCount, out fractionalParts);
            someParents.CopyTo(allParents, 0);

            // Remainder phase
            int extraParentsCount = totalParentsCount - someParents.Count;
            Individual<T>[] extraParents = RemainderPhase(population, fractionalParts, extraParentsCount);
            Array.Copy(extraParents, 0, allParents, someParents.Count, extraParentsCount);

            // Shuffle the parents to avoid identical ones being grouped together
            RandomNumberGenerationUtilities.Shuffle<Individual<T>>(allParents);

            // Place the parents in groups
            int counter = 0;
            var parentGroups = new Individual<T>[parentGroupsCount][];
            for (int group = 0; group < parentGroupsCount; ++group)
            {
                parentGroups[group] = new Individual<T>[parentsPerGroup];
                for (int parent = 0; parent < parentsPerGroup; ++parent)
                {
                    parentGroups[group][parent] = allParents[counter];
                    ++counter;
                }
            }
            return parentGroups;
        }

        // For each Individual: expectation = totalParentsCount * expectation / sum(expectations)
        private LinkedList<Individual<T>> IntegralPhase(Individual<T>[] population, int totalParentsCount,
            out double[] fractionalParts)
        {
            // Normalize expectations so that they are positive and their sum = totalParentsCount
            double[] expectations = fitnessScaling.CalculateExpectations(population);
            double scale = totalParentsCount / expectations.Sum();
            for (int i = 0; i < expectations.Length; ++i)
            {
                expectations[i] *= scale;
            }

            var selectedParents = new LinkedList<Individual<T>>();
            fractionalParts = new double[population.Length];
            for (int parent = 0; parent < population.Length; ++parent)
            {
                // Each parent will be selected times equal to the integer part of its scaled expectation.
                int integral = (int)expectations[parent];
                for (int i = 0; i < integral; ++i) selectedParents.AddLast(population[parent]);

                // The fractional part will be used later
                fractionalParts[parent] = expectations[parent] - integral;
            }

            return selectedParents;
        }

        /// <summary>
        /// Template method.
        /// </summary>
        /// <param name="population"></param>
        /// <param name="fractionalParts"></param>
        /// <param name="extrasParentsCount"></param>
        /// <returns></returns>
        protected abstract Individual<T>[] RemainderPhase(Individual<T>[] population, double[] fractionalParts,
            int extrasParentsCount);
    }
}
