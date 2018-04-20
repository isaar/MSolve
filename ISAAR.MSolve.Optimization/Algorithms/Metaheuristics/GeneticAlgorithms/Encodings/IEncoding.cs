using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings
{
    public interface IEncoding<T>
    {
        T[] ComputeGenotype(double[] phenotype);
        double[] ComputePhenotype(T[] genotype);
        //int[] IntegerPhenotype(T[] genotype);
    }
}
