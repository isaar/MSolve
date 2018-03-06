using ISAAR.MSolve.Optimization.Commons;
using ISAAR.MSolve.Optimization.Problems;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations
{
    public class UniformMutation: IMutationStrategy<double>
    {
        private readonly double mutationProbability;
        private readonly IGenerator rng;
        private readonly double[] lowerBounds;
        private readonly double[] upperBounds;

        public UniformMutation(OptimizationProblem problem, double mutationProbability) : 
                            this(problem, mutationProbability, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public UniformMutation(OptimizationProblem problem, double mutationProbability, IGenerator randomNumberGenerator)
        {
            //problem.CheckInput();
            this.lowerBounds = problem.LowerBound;
            this.upperBounds = problem.UpperBound;

            if (mutationProbability < 0 || mutationProbability > 1)
            {
                throw new ArgumentException("The mutation probability of each gene must belong to the interval [0,1], but was "
                                             + mutationProbability);
            }
            this.mutationProbability = mutationProbability;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public void Apply(Individual<double>[] population)
        {
            CanonicalVersion(population);
            //FastVersion(population);
        }

        // Running time = O(populationSize * chromosomeSize). Not terrible for real or integer value encodings
        private void CanonicalVersion(Individual<double>[] population)
        {
            int genesCount = population[0].Chromosome.Length; //No checking for the other chromosomes
            foreach (var individual in population)
            {
                double[] chromosome = individual.Chromosome;
                for (int gene = 0; gene < genesCount; ++gene)
                {
                    if (rng.NextDouble() < mutationProbability)
                    {
                        chromosome[gene] = rng.NextDouble(lowerBounds[gene], upperBounds[gene]);
                    }
                }
            }
        }

        // Faster than canonical version, but I am not sure if it is valid statistically. 
        // The canonical version performs (populationSize * chromosomeSize * mutationProbability) mutations on average. 
        // In contrast this version always performs exactly as many mutations. 
        // Also it is possible for 2 or more mutations to occur on the same gene of the same chromosome, thus negating each other.
        private void FastVersion(Individual<double>[] population)
        {
            int genesCount = population[0].Chromosome.Length; //No checking for the other chromosomes
            int totalMutations = (int)(population.LongLength * (mutationProbability * genesCount));
            for (int i = 0; i < totalMutations; ++i)
            {
                var chromosome = population[rng.Next(population.Length)].Chromosome;
                int gene = rng.Next(genesCount);
                chromosome[gene] = rng.NextDouble(lowerBounds[gene], upperBounds[gene]);
            }
        }
    }
}
