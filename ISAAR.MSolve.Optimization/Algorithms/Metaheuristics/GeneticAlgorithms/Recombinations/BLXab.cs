using ISAAR.MSolve.Optimization.Commons;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations
{
    public class BLXab : IRecombinationStrategy<double>
    {
        private readonly IGenerator rng;
        private readonly double alpha;
        private readonly double beta;

        public BLXab(double alpha = 0.75, double beta = 0.25) :
            this(alpha, beta, RandomNumberGenerationUtilities.troschuetzRandom)
        { }

        public BLXab(double alpha, double beta, IGenerator randomNumberGenerator)
        {
            // TODO: find min, max for alpha and beta
            if (alpha <= beta)
            {
                throw new ArgumentException("Parameter alpha = " + alpha + " must be greater than parameter beta = " + beta);
            }
            this.alpha = alpha;
            this.beta = beta;

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
                if (parent0.Fitness <= parent1.Fitness)
                {
                    offsprings[i] = new Individual<double>(Blend(parent0.Chromosome, parent1.Chromosome));
                }
                else
                {
                    offsprings[i] = new Individual<double>(Blend(parent1.Chromosome, parent0.Chromosome));
                }
            }
            return offsprings;
        }

        private double[] Blend(double[] bestParent, double[] worstParent)
        {
            int genesCount = bestParent.Length; //No checking for the other chromosomes
            double[] offspring = new double[genesCount];
            for (int gene = 0; gene < genesCount; ++gene)
            {
                double range = Math.Abs(bestParent[gene] - worstParent[gene]);
                if (bestParent[gene] <= worstParent[gene])
                {
                    double min = bestParent[gene] - alpha * range;
                    double max = worstParent[gene] + beta * range;
                    offspring[gene] = rng.NextDouble(min, max);
                }
                else
                {
                    double min = worstParent[gene] - beta * range;
                    double max = bestParent[gene] + alpha * range;
                    offspring[gene] = rng.NextDouble(min, max);
                }
            }
            return offspring;
        }
    }
}
