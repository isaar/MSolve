using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Initialization.HaltonPoints
{
    // It is not thread safe yet!
    public class HaltonPointGenerator
    {
        private readonly int dimension;
        private readonly HaltonSequence[] sequences;
        private int currentIndex;

        public HaltonPointGenerator(int dimension)
        {
            if (dimension < 1)
            {
                throw new ArgumentException("Dimension must be >= 1, but was " + dimension);
            }
            this.dimension = dimension;
            this.sequences = new HaltonSequence[dimension];

            int[] primes = SelectPrimeBases(dimension);
            for (int i = 0; i < dimension; ++i)
            {
                sequences[i] = new HaltonSequence(primes[i]);
            }

            currentIndex = 1;
        }

        public double[] NextPoint()
        {
            double[] point = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                point[i] = sequences[i].ElementAt(currentIndex);
            }
            ++currentIndex;
            return point;
        }

        private int[] SelectPrimeBases(int dimension)
        {
            IPrimeGenerator primeGenerator = new First50PrimesGenerator();
            int[] firstPrimes = primeGenerator.FirstPrimes(dimension);
            return firstPrimes; // there may be a correlation for high dimensions
        }
    }
}
