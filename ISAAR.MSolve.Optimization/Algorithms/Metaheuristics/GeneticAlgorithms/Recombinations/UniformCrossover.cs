using ISAAR.MSolve.Optimization.Commons;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Troschuetz.Random;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations
{
    public class UniformCrossover<T> : IRecombinationStrategy<T>
    {
        private readonly IGenerator rng;

        public UniformCrossover() : this(RandomNumberGenerationUtilities.troschuetzRandom) { }

        public UniformCrossover(IGenerator randomNumberGenerator)
        {
            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public Individual<T>[] Apply(ISelectionStrategy<T> selection, Individual<T>[] population, int offspringsCount)
        {
            // Select the necessary parents
            int pairsCount = (offspringsCount - 1) / 2 + 1;
            var parentGroups = selection.Apply(population, pairsCount, 2);

            // Create the offsprings
            var offsprings = new Individual<T>[offspringsCount];
            for (int pair = 0; pair < pairsCount; ++pair)
            {
                T[] offspring1, offspring2;
                Crossover(parentGroups[pair][0].Chromosome, parentGroups[pair][1].Chromosome, out offspring1, out offspring2);
                offsprings[2 * pair] = new Individual<T>(offspring1);
                if (2 * pair + 1 < offspringsCount) offsprings[2 * pair + 1] = new Individual<T>(offspring2);
            }
            return offsprings;
        }

        private void Crossover(T[] parent1, T[] parent2, out T[] offspring1, out T[] offspring2)
        {
            int genesCount = parent1.Length; //No checking for the other chromosomes
            offspring1 = new T[genesCount];
            offspring2 = new T[genesCount];
            for (int gene = 0; gene < genesCount; ++gene)
            {
                if (rng.NextBoolean()) // true => crossover
                {
                    offspring1[gene] = parent2[gene];
                    offspring2[gene] = parent1[gene];
                }
                else // false => inherit from corresponding parent
                {
                    offspring1[gene] = parent1[gene];
                    offspring2[gene] = parent2[gene];
                }
            }
        }
    }
}
