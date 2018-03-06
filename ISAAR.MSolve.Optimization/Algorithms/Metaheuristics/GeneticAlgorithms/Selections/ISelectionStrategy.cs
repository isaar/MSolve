using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections
{
    public interface ISelectionStrategy<T>
    {
        Individual<T>[][] Apply(Individual<T>[] population, int parentGroupsCount, int parentsPerGroup);
    }
}
