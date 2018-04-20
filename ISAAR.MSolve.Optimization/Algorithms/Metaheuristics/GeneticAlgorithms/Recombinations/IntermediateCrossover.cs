using ISAAR.MSolve.Optimization.Commons;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations
{
    public class IntermediateCrossover : IRecombinationStrategy<double>
    {
        private readonly IGenerator rng;
        private readonly double firstParentWeight;

        public IntermediateCrossover(double firstParentWeight = 1.0) : 
            this(firstParentWeight, RandomNumberGenerationUtilities.troschuetzRandom) { }

        public IntermediateCrossover(double firstParentWeight, IGenerator randomNumberGenerator)
        {
            if ((firstParentWeight < 0) || (firstParentWeight > 1))
            {
                throw new ArgumentException("The weight coefficient of the first parent selected for recombination"
                    + " must belong to the interval [0, 1], but was " + firstParentWeight);
            }
            this.firstParentWeight = firstParentWeight;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public Individual<double>[] Apply(ISelectionStrategy<double> selection, Individual<double>[] population, int offspringsCount)
        {
            // Select the necessary parents
            int pairsCount = (offspringsCount - 1) / 2 + 1;
            var parentGroups = selection.Apply(population, pairsCount, 2);

            // Create the offsprings
            var offsprings = new Individual<double>[offspringsCount];
            for (int pair = 0; pair < pairsCount; ++pair)
            {
                double[] offspring1, offspring2;
                Interpolate(parentGroups[pair][0].Chromosome, parentGroups[pair][1].Chromosome, out offspring1, out offspring2);
                offsprings[2 * pair] = new Individual<double>(offspring1);
                if (2 * pair + 1 < offspringsCount) offsprings[2 * pair + 1] = new Individual<double>(offspring2);
            }
            return offsprings;
        }

        private void Interpolate(double[] parent1, double[] parent2, out double[] offspring1, out double[] offspring2)
        {
            int genesCount = parent1.Length; //No checking for the other chromosomes
            offspring1 = new double[genesCount];
            offspring2 = new double[genesCount];
            double weight = rng.NextDouble() * firstParentWeight;
            for (int gene = 0; gene < genesCount; ++gene)
            {
                offspring1[gene] = weight * parent1[gene] + (1 - weight) * parent2[gene];
                offspring2[gene] = weight * parent2[gene] + (1 - weight) * parent1[gene];
            }
        }
    }
}
