using ISAAR.MSolve.Optimization.Commons;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations
{
    public class BitFlipMutation : IMutationStrategy<bool>
    {
        private readonly double mutationProbability;
        private readonly IGenerator rng;

        public BitFlipMutation(double mutationProbability) : 
                            this(mutationProbability, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public BitFlipMutation(double mutationProbability, IGenerator randomNumberGenerator)
        {
            if (mutationProbability < 0 || mutationProbability > 1)
            {
                throw new ArgumentException("The mutation probability of each gene must belong to the interval [0,1], but was "
                                             + mutationProbability);
            }
            this.mutationProbability = mutationProbability;
            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public void Apply(Individual<bool>[] population)
        {
            //CanonicalVersion(population);
            FastVersion(population);
        }

        // Running time = O(populationSize * chromosomeSize). Particularly slow for binary encoding, where chromosomeSize is large.
        private void CanonicalVersion(Individual<bool>[] population)
        {
            int genesCount = population[0].Chromosome.Length; //No checking for the other chromosomes
            foreach (var individual in population)
            {
                bool[] chromosome = individual.Chromosome;
                for (int gene = 0; gene < genesCount; ++gene)
                {
                    if (rng.NextDouble() < mutationProbability) chromosome[gene] = !chromosome[gene];
                }
            }
        }

        // Faster than canonical version, but I am not sure if it is valid statistically. 
        // The canonical version performs (populationSize * chromosomeSize * mutationProbability) mutations on average. 
        // In contrast this version always performs exactly as many mutations. 
        // Also it is possible for 2 or more mutations to occur on the same gene of the same chromosome, thus negating each other.
        private void FastVersion(Individual<bool>[] population)
        {
            int genesCount = population[0].Chromosome.Length; //No checking for the other chromosomes
            int totalMutations = (int)(population.LongLength * (mutationProbability * genesCount));
            for (int i = 0; i < totalMutations; ++i)
            {
                var chromosome = population[rng.Next(population.Length)].Chromosome;
                int gene = rng.Next(genesCount);
                chromosome[gene] = !chromosome[gene];
            }
        }
    }
}
