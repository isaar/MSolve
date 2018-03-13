using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.PopulationStrategies
{
    public interface IPopulationStrategy<T>
    {
        Individual<T>[] CreateNextGeneration(Individual<T>[] originalPopulation, ISelectionStrategy<T> selection, 
                                      IRecombinationStrategy<T> recombination, IMutationStrategy<T> mutation); 
    }
}
