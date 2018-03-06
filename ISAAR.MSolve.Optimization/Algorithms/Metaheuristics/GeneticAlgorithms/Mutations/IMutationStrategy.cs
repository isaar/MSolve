using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations
{
    public interface IMutationStrategy<T>
    {
        void Apply(Individual<T>[] population);
    }
}
