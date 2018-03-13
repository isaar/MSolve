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
    public class UnfitEliminationStrategy<T> : IPopulationStrategy<T>
    {
        private readonly int populationSize;
        private readonly int elitesCount;
        private readonly int offspringsCount;
        private readonly int survivorsCount; //excluding Elites
        private readonly IGenerator rng;

        // Need to find a better way to check elitesCount against populationSize
        public UnfitEliminationStrategy(int populationSize, int elitesCount, double eliminationPercentage) : 
                      this(populationSize, elitesCount, eliminationPercentage, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public UnfitEliminationStrategy(int populationSize, int elitesCount, double eliminationPercentage, 
                                        IGenerator randomNumberGenerator)
        {
            if (populationSize < 1) throw new ArgumentException("There population size must be >= 1");
            this.populationSize = populationSize;

            if ((elitesCount < 0) || (elitesCount >= populationSize))
            {
                throw new ArgumentException("The number of elites must belong to the interval [0, populationSize-1), but was "
                    + elitesCount);
            }
            this.elitesCount = elitesCount;

            int eliminationsCount = (int)Math.Round(eliminationPercentage * populationSize);
            if ((eliminationsCount < 0) || (eliminationsCount >= populationSize))
            {
                throw new ArgumentException("The number of individuals eliminated before selection must belong to the interval" 
                    + " [0, populationSize-1), but was " + elitesCount);
            }
            this.offspringsCount = eliminationsCount;
            this.survivorsCount = populationSize - eliminationsCount - elitesCount;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public Individual<T>[] CreateNextGeneration(Individual<T>[] originalPopulation, ISelectionStrategy<T> selection,
                                                    IRecombinationStrategy<T> recombination, IMutationStrategy<T> mutation)
        {
            Array.Sort(originalPopulation);

            // Selection and rcombination: Only the elites and the rest of survivors will have a chance to reproduce
            Individual<T>[] selectionPool = new Individual<T>[elitesCount + survivorsCount];
            Array.Copy(originalPopulation, selectionPool, selectionPool.Length);
            // TODO: Redundant copying. A linked list would be better.
            Individual<T>[] offsprings = recombination.Apply(selection, selectionPool, offspringsCount);

            // Mutation will be applied to the survivors and their offsprings but not on the elites
            Individual<T>[] mutants = new Individual<T>[survivorsCount + offsprings.Length];
            Array.Copy(originalPopulation, elitesCount, mutants, 0, survivorsCount);
            Array.Copy(offsprings, 0, mutants, survivorsCount, offspringsCount);
            mutation.Apply(mutants);

            // The next population will contain the elites, the mutated survivors and the mutated offsprings
            Individual<T>[] nextPopulation = new Individual<T>[originalPopulation.Length];
            Array.Copy(originalPopulation, nextPopulation, elitesCount);
            Array.Copy(mutants, 0, nextPopulation, elitesCount, survivorsCount + offspringsCount);
            return nextPopulation;
        }
    }
}
