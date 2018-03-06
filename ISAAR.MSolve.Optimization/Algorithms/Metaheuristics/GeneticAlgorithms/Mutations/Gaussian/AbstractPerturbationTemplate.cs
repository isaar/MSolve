using ISAAR.MSolve.Optimization.Commons;
using ISAAR.MSolve.Optimization.Problems;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations.Gaussian
{
    public abstract class AbstractPerturbationTemplate : IMutationStrategy<double>
    {
        public void Apply(Individual<double>[] population)
        {
            foreach (var individual in population)
            {
                double[] chromosome = individual.Chromosome;
                individual.Chromosome = VectorOperations.Add(individual.Chromosome, ComputePerturbations());
            }
        }

        protected abstract double[] ComputePerturbations();
    }
}
