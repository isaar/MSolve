using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings
{
    class RealCoding : IEncoding<double>
    {
        public double[] ComputeGenotype(double[] phenotype)
        {
            return phenotype;
        }

        public double[] ComputePhenotype(double[] genotype)
        {
            return genotype;
        }
    }
}
