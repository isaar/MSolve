using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.PopulationStrategies;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;
using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Initialization;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Problems;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms
{
    public partial class GeneticAlgorithm<T> : IOptimizationAlgorithm
    {
        public abstract class Builder
        {
            #region fields, properties, constructor
            // Optim problem fields.
            OptimizationProblem problem;

            protected Builder(OptimizationProblem problem)
            {
                this.problem = problem;
            }

            // General optimiazation algorithm parameterss
            public IOptimizationLogger Logger { get; set; }
            public IConvergenceCriterion ConvergenceCriterion { get; set; }
            public IInitializer<double> Initializer { get; set; }

            // GA parameters
            public IEncoding<T> Encoding { get; set; }
            public int PopulationSize { get; set; }
            public IPopulationStrategy<T> PopulationStrategy { get; set; }
            public ISelectionStrategy<T> Selection { get; set; }
            public IRecombinationStrategy<T> Recombination { get; set; }
            public IMutationStrategy<T> Mutation { get; set; }

            // Default GA parameters -> msut be implemented
            protected abstract IEncoding<T> DefaultEncoding { get; }
            protected abstract int DefaultPopulationSize { get; }
            protected abstract IPopulationStrategy<T> DefaultPopulationStrategy { get; }
            protected abstract ISelectionStrategy<T> DefaultSelection { get; }
            protected abstract IRecombinationStrategy<T> DefaultRecombination { get; }
            protected abstract IMutationStrategy<T> DefaultMutation { get; }

            #endregion

            #region methods
            public GeneticAlgorithm<T> BuildAlgorithm()
            {

                CheckUserParameters();
                ApplyDefaultParameters();
                return new GeneticAlgorithm<T>(problem.Dimension, 0, problem.DesignFactory, 
                                               Logger, ConvergenceCriterion, Initializer, 
                                               Encoding, PopulationSize, PopulationStrategy, Selection, Recombination, Mutation);
            }

            private void ApplyDefaultParameters()
            {
                if (Logger == null) Logger = new BestOfIterationLogger();
                if (ConvergenceCriterion == null) ConvergenceCriterion = new MaxIterations(100 * problem.Dimension);
                if (Initializer == null) Initializer = new RealUniformRandomInitializer(problem);

                // The GA specific default parameters are provided be the implementation of this Builder
                if (Encoding == null) Encoding = DefaultEncoding;
                if (PopulationSize == 0) PopulationSize = DefaultPopulationSize;
                if (PopulationStrategy == null) PopulationStrategy = DefaultPopulationStrategy;
                if (Selection == null) Selection = DefaultSelection;
                if (Recombination == null) Recombination = DefaultRecombination;
                if (Mutation == null) Mutation = DefaultMutation;
            }

            // Crash when user provides incompatible parameters. Warn him when the parameters, although legal, may result in poor performance 
            // Will ignore the default values. TODO find a better way to handle this.
            private void CheckUserParameters()
            {
                if ((PopulationSize != 0) && (PopulationSize < 1))
                {
                    throw new ArgumentException("Population size must be at least 1, but was " + PopulationSize);
                }
            }
            #endregion
        }
    }
}
