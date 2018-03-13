using ISAAR.MSolve.Optimization.Commons;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations
{
    
    public class BLXa: IRecombinationStrategy<double>
    {
        private readonly IGenerator rng;
        private readonly double alpha;

        public BLXa(double alpha = 0.5) : 
            this(alpha, RandomNumberGenerationUtilities.troschuetzRandom) { }

        public BLXa(double alpha, IGenerator randomNumberGenerator)
        {
            // TODO: find min, max alpha
            //if (true)
            //{
            //    throw new ArgumentException("The alpha parameter must belong to the interval [], but was " + alpha);
            //}
            this.alpha = alpha;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public Individual<double>[] Apply(ISelectionStrategy<double> selection, Individual<double>[] population, int offspringsCount)
        {
            // Select the necessary parents
            var parentGroups = selection.Apply(population, offspringsCount, 2);

            // Create the offsprings
            var offsprings = new Individual<double>[offspringsCount];
            for (int i = 0; i < offspringsCount; ++i)
            {
                var parent0 = parentGroups[i][0];
                var parent1 = parentGroups[i][1];
                offsprings[i] = new Individual<double>(Blend(parent0.Chromosome, parent1.Chromosome));
            }
            return offsprings;
        }

        private double[] Blend(double[] parent1, double[] parent2)
        {
            int genesCount = parent1.Length; // No checking for the other chromosomes
            double[] offspring = new double[genesCount];
            for (int gene = 0; gene < genesCount; ++gene)
            {
                double range = Math.Abs(parent1[gene] - parent2[gene]);
                double min = Math.Min(parent1[gene], parent2[gene]) - alpha * range;
                double max = Math.Max(parent1[gene], parent2[gene]) + alpha * range;
                offspring[gene] = rng.NextDouble(min, max);
            }
            return offspring;
        }
    }
}
