using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.PopulationStrategies
{
    class SteadyStateStrategy<T> : IPopulationStrategy<T>
    {
        public Individual<T>[] CreateNextGeneration(Individual<T>[] originalPopulation, ISelectionStrategy<T> selection, 
                                                    IRecombinationStrategy<T> recombination, IMutationStrategy<T> mutation)
        {
            throw new NotImplementedException();
        }
    }
}
