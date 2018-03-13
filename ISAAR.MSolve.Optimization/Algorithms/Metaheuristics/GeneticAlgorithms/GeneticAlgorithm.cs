using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.PopulationStrategies;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;
using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Problems;

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Optimization.Initialization;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms
{
    /// <summary>
    /// Genetic Algorithm with non adaptive parameters
    /// </summary>
    public partial class GeneticAlgorithm<T> : IOptimizationAlgorithm
    {
        #region fields, properties, constructor
        // Optim problem fields
        private readonly int continuousVariablesCount;
        private readonly int integerVariablesCount;
        private readonly IDesignFactory designFactory;

        // General optim algorithm params
        private readonly IOptimizationLogger logger;
        private readonly IConvergenceCriterion convergenceCriterion;
        private readonly IInitializer<double> initializer;

        // GA params
        private readonly IEncoding<T> encoding;
        private readonly int populationSize;
        private readonly IPopulationStrategy<T> populationStrategy;
        private readonly ISelectionStrategy<T> selection;
        private readonly IRecombinationStrategy<T> recombination;
        private readonly IMutationStrategy<T> mutation;

        private Individual<T>[] population;

        private GeneticAlgorithm(int continuousVariablesCount, int integerVariablesCount, IDesignFactory designFactory,
            IOptimizationLogger logger, IConvergenceCriterion convergenceCriterion, IInitializer<double> initializer,
            IEncoding<T> encoding, int populationSize, IPopulationStrategy<T> populationStrategy,
            ISelectionStrategy<T> selection, IRecombinationStrategy<T> recombination, IMutationStrategy<T> mutation)
        {
            this.continuousVariablesCount = continuousVariablesCount;
            this.integerVariablesCount = integerVariablesCount;
            this.designFactory = designFactory;

            this.logger = logger;
            this.convergenceCriterion = convergenceCriterion;
            this.initializer = initializer;

            this.encoding = encoding;
            this.populationSize = populationSize;
            this.populationStrategy = populationStrategy;
            this.selection = selection;
            this.recombination = recombination;
            this.mutation = mutation;

            this.CurrentIteration = -1; // Initialization phase is not counted towards iterations
            this.BestPosition = null;
            this.BestFitness = double.MaxValue;
        }

        public double BestFitness { get; private set; }
        public double[] BestPosition { get; private set; }
        public int CurrentIteration { get; private set; }
        public double CurrentFunctionEvaluations { get; private set; }
        #endregion

        #region public methods
        public void Solve()
        {
            CurrentIteration = 0;
            Initialize();
            logger.Log(this);

            while (!convergenceCriterion.HasConverged(this))
            {
                ++CurrentIteration;
                Iterate();
                logger.Log(this);
            }
        }
        #endregion

        #region private methods
        private void Initialize()
        {
            population = new Individual<T>[populationSize];
            for (int i = 0; i < populationSize; ++i)
            {
                double[] phenotype = initializer.Generate();
                population[i] = new Individual<T>(encoding.ComputeGenotype(phenotype));
            }

            EvaluateCurrentIndividuals();
        }

        private void Iterate()
        {
            population = populationStrategy.CreateNextGeneration(population, selection, recombination, mutation);

            EvaluateCurrentIndividuals();
        }

        private void EvaluateCurrentIndividuals()
        {
            foreach (Individual<T> individual in population)
            {
                if (!individual.IsEvaluated)
                {
                    double[] phenotype = encoding.ComputePhenotype(individual.Chromosome);
                    IDesign design = designFactory.CreateDesign(phenotype);
                    individual.Fitness = design.ObjectiveValues[0];
                }
            }
            CurrentFunctionEvaluations += populationSize;
            UpdateBest();
        }

        private void UpdateBest()
        {
            foreach (Individual<T> individual in population)
            {
                if (individual.Fitness < BestFitness)
                {
                    BestFitness = individual.Fitness;
                    BestPosition = encoding.ComputePhenotype(individual.Chromosome); // Redundant conversion; must be removed 
                }
            }
        }
        #endregion
    }
}