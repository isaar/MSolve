using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations
{
    public interface IRecombinationStrategy<T>
    {
        Individual<T>[] Apply(ISelectionStrategy<T> selection, Individual<T>[] population, int offspringsCount);
    }
}
