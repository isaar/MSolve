using ISAAR.MSolve.Optimization.Commons;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections
{
    class UniformRandomSelection<T>: ISelectionStrategy<T>
    {
        private readonly bool allowIdenticalParents;
        private readonly IGenerator rng;

        public UniformRandomSelection(bool allowIdenticalParents = false) : 
            this(allowIdenticalParents, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public UniformRandomSelection(bool allowIdenticalParents, IGenerator randomNumberGenerator)
        {
            this.allowIdenticalParents = allowIdenticalParents;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public Individual<T>[][] Apply(Individual<T>[] population, int parentGroupsCount, int parentsPerGroup)
        {
            var parentGroups = new Individual<T>[parentGroupsCount][];
            for (int group = 0; group < parentGroupsCount; ++group)
            {
                parentGroups[group] = new Individual<T>[parentsPerGroup];
                for (int parent = 0; parent < parentsPerGroup; ++parent)
                {
                    Individual<T> individual = population[rng.Next(population.Length)];
                    if (!allowIdenticalParents)
                    {
                        while (parentGroups[group].Contains<Individual<T>>(individual))
                        {
                            individual = population[rng.Next(population.Length)];
                        }
                    }
                    parentGroups[group][parent] = individual;
                }
            }
            return parentGroups;
        }
    }
}
