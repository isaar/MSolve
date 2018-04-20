using ISAAR.MSolve.Optimization.Commons;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Troschuetz.Random;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections
{
    public class RouletteWheelSelection<T> : ISelectionStrategy<T>
    {
        private readonly IFitnessScalingStrategy<T> fitnessScaling;
        private readonly bool allowIdenticalParents;
        private readonly IGenerator rng;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="fitnessScaling"></param>
        /// <param name="allowIdenticalParents">If set to false and the fitness scaling is too biased towards the fittest 
        ///     individuals, performance will degrade dramatically due to excessive sampling</param>
        public RouletteWheelSelection(IFitnessScalingStrategy<T> fitnessScaling, bool allowIdenticalParents = true) :
            this(fitnessScaling, allowIdenticalParents, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public RouletteWheelSelection(IFitnessScalingStrategy<T> fitnessScaling, bool allowIdenticalParents, 
            IGenerator randomNumberGenerator)
        {
            if (fitnessScaling == null) throw new ArgumentException("The fitness scaling strategy must not be null");
            this.fitnessScaling = fitnessScaling;

            this.allowIdenticalParents = allowIdenticalParents;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public Individual<T>[][] Apply(Individual<T>[] population, int parentGroupsCount, int parentsPerGroup)
        {
            double[] expectations = fitnessScaling.CalculateExpectations(population);
            Roulette roulette = Roulette.CreateFromPositive(expectations, rng);

            var parentGroups = new Individual<T>[parentGroupsCount][];
            for (int group = 0; group < parentGroupsCount; ++group)
            {
                parentGroups[group] = new Individual<T>[parentsPerGroup];
                for (int parent = 0; parent < parentsPerGroup; ++parent)
                {
                    Individual<T> individual = population[roulette.SpinWheelWithBall()];
                    if (!allowIdenticalParents)
                    {
                        while (parentGroups[group].Contains<Individual<T>>(individual))
                        {
                            individual = population[roulette.SpinWheelWithBall()];
                        }
                    }
                    parentGroups[group][parent] = individual;
                }
            }
            return parentGroups;
        }
    }
}
