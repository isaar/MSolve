using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling;
using ISAAR.MSolve.Optimization.Commons;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections
{
    public class StochasticUniversalSampling<T> : ISelectionStrategy<T>
    {
        private readonly IFitnessScalingStrategy<T> fitnessScaling;
        private readonly IGenerator rng;

        public StochasticUniversalSampling(IFitnessScalingStrategy<T> fitnessScaling) :
            this(fitnessScaling, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public StochasticUniversalSampling(IFitnessScalingStrategy<T> fitnessScaling, IGenerator randomNumberGenerator)
        {
            if (fitnessScaling == null) throw new ArgumentException("The fitness scaling strategy must not be null");
            this.fitnessScaling = fitnessScaling;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public Individual<T>[][] Apply(Individual<T>[] population, int parentGroupsCount, int parentsPerGroup)
        {
            Individual<T>[] allParents = SampleAllAtOnce(population, parentGroupsCount * parentsPerGroup);
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

        private Individual<T>[] SampleAllAtOnce(Individual<T>[] population, int totalParentsCount)
        {
            double[] expectations = fitnessScaling.CalculateExpectations(population);
            Roulette roulette = Roulette.CreateFromPositive(expectations, rng);
            int[] selectedIndexes = roulette.SpinWheelWithPointers(totalParentsCount);

            var allParents = new Individual<T>[totalParentsCount];
            for (int i = 0; i < totalParentsCount; ++i)
            {
                allParents[i] = population[selectedIndexes[i]];
            }
            return allParents;
        }
    }
}
