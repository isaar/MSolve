using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Optimization.Commons;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections
{
    public class TournamentSelection<T> : ISelectionStrategy<T>
    {
        private readonly int tournamentSize;
        private readonly bool sampleWithReplacement;
        private readonly bool allowIdenticalParents;
        private readonly IGenerator rng;

        public TournamentSelection(int tournamentSize, bool sampleWithReplacement = true, bool allowIdenticalParents = false):
            this(tournamentSize, sampleWithReplacement, allowIdenticalParents, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tournamentSize">The number of individuals in each tournament. It must be greater than 1 and 
        ///     less or equal to the population size. For tournament size = 1 use <see cref="UniformRandomSelection{T}"/> 
        ///     instead.</param>
        /// <param name="sampleWithReplacement">If set to false, winners of tournaments will be removed from the selection pool 
        ///     for that generation, meaning that each individual will be selected once at most. Caution: If set to false, 
        ///     the total number of parents cannot exceed the population size and the recombination strategy cannot request more 
        ///     parents than that. In this case there won't be identical parents in the same parent group, thus opting to allow 
        ///     identical parents will have no effect.</param>
        /// <param name="allowIdenticalParents">If set to false, sampling will be repeated until all parents of the same group are 
        ///     uniqe. This may degrade performance, especially with large tournament sizes.</param>
        /// <param name="randomNumberGenerator">A random number generator.</param>
        public TournamentSelection(int tournamentSize, bool sampleWithReplacement, bool allowIdenticalParents, 
            IGenerator randomNumberGenerator)
        {
            if (tournamentSize < 2)
            {
                throw new ArgumentException("The tournament size must be at least 2, but was " + tournamentSize);
            }
            this.tournamentSize = tournamentSize;

            this.sampleWithReplacement = sampleWithReplacement;
            this.allowIdenticalParents = allowIdenticalParents;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public Individual<T>[][] Apply(Individual<T>[] population, int parentGroupsCount, int parentsPerGroup)
        {
            if (tournamentSize > population.Length)
            {
                throw new ArgumentException("The tournament size: " + tournamentSize + 
                    " must be at lower than or equal to the population size: " + population.Length);
            }
            if (!sampleWithReplacement)
            {
                if (parentGroupsCount * parentsPerGroup > population.Length)
                {
                    throw new ArgumentException("Tournament selection cannot generate " + parentGroupsCount + " * "
                        + parentsPerGroup + " parents without replacement");
                }
            }

            var selectionPool = new Bag<Individual<T>>(population);
            var parentGroups = new Individual<T>[parentGroupsCount][];
            for (int group = 0; group < parentGroupsCount; ++group)
            {
                parentGroups[group] = new Individual<T>[parentsPerGroup];
                for (int parent = 0; parent < parentsPerGroup; ++parent)
                {
                    Individual<T> individual = ConductTournament(selectionPool);
                    if (!allowIdenticalParents)
                    {
                        while (parentGroups[group].Contains<Individual<T>>(individual))
                        {
                            individual = ConductTournament(selectionPool);
                        }
                    }
                    parentGroups[group][parent] = individual;
                }
            }
            return parentGroups;
        }

        private Individual<T> ConductTournament(Bag<Individual<T>> selectionPool)
        {
            Individual<T> winner = selectionPool.GetRandom();
            for (int i = 1; i < tournamentSize; ++i)
            {
                Individual<T> competitor = selectionPool.GetRandom();
                if (competitor.Fitness < winner.Fitness) winner = competitor;
            }
            if (!sampleWithReplacement) selectionPool.Remove(winner);
            return winner;
        }
    }
}
