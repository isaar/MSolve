using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.PopulationStrategies;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms
{
    public class BinaryGeneticAlgorithmBuilder : GeneticAlgorithm<bool>.Builder
    {
        private OptimizationProblem problem;

        public BinaryGeneticAlgorithmBuilder(OptimizationProblem problem) : base(problem)
        {
            this.problem = problem;
        }

        protected override IEncoding<bool> DefaultEncoding
        {
            get // good for continuous design variables 
            {
                return new GrayCoding(problem, 32, 8); // sizes of float, char 
            }
        }

        protected override IMutationStrategy<bool> DefaultMutation
        {
            get // arbitrary
            {
                return new BitFlipMutation(0.2);
            }
        }

        protected override int DefaultPopulationSize
        {
            get // Matlab defaults
            {
                return (problem.Dimension <= 5) ? 50 : 200;

                //if (integerVariablesCount == 0) // continuous problem
                //{
                //    PopulationSize = (continuousVariablesCount <= 5) ? 50 : 200;
                //}
                //else // Mixed integer problem
                //{
                //    PopulationSize = Math.Min(100, Math.Max(40, 10 * (continuousVariablesCount + integerVariablesCount)));
                //}
            }
        }

        protected override IPopulationStrategy<bool> DefaultPopulationStrategy
        {
            get // arbitrary strategy with Matlab's default elitism 
            {
                return new StandardPopulationStrategy<bool>(PopulationSize, (int)Math.Round(0.05 * PopulationSize));
            }
        }

        protected override IRecombinationStrategy<bool> DefaultRecombination
        {
            get // Matlab defaults
            {
                return new UniformCrossover<bool>();
            }
        }

        protected override ISelectionStrategy<bool> DefaultSelection
        {
            get // Matlab defaults
            {
                return new RouletteWheelSelection<bool>(new InverseRankScaling<bool>(0.5));
            }
        }
    }
}
