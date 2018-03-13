using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;
using Troschuetz.Random;
using ISAAR.MSolve.Optimization.Commons;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.PopulationStrategies
{
    class StandardPopulationStrategy<T> : IPopulationStrategy<T>
    {
        private readonly int populationSize;
        private readonly int elitesCount;
        private readonly IGenerator rng;
        
        // Need to find a better way to check elitesCount against populationSize
        public StandardPopulationStrategy(int populationSize, int elitesCount) : 
                            this(populationSize, elitesCount, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public StandardPopulationStrategy(int populationSize, int elitesCount, IGenerator randomNumberGenerator)
        {
            if (populationSize < 1) throw new ArgumentException("There population size must be >= 1");
            this.populationSize = populationSize;

            if ((elitesCount < 0) || (elitesCount >= populationSize))
            {
                throw new ArgumentException("The number of elites must belong to the interval [0, populationSize-1), but was "
                    + elitesCount);
            }
            this.elitesCount = elitesCount;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public Individual<T>[] CreateNextGeneration(Individual<T>[] originalPopulation, ISelectionStrategy<T> selection, 
                                             IRecombinationStrategy<T> recombination, IMutationStrategy<T> mutation)
        {
            Array.Sort(originalPopulation); // sorting may not always be mandatory (e.g. 1-3 elites, selection does not sort)
            int offspringsCount = populationSize - elitesCount;
            // TODO: Redundant copying. A linked list would be better.
            Individual<T>[] offsprings = recombination.Apply(selection, originalPopulation, offspringsCount);
            mutation.Apply(offsprings);

            Individual<T>[] nextPopulation = new Individual<T>[populationSize];
            Array.Copy(originalPopulation, nextPopulation, elitesCount);
            Array.Copy(offsprings, 0, nextPopulation, elitesCount, offsprings.Length);
            return nextPopulation;
        }
    }
}
