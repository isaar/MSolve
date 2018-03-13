using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Initialization.HaltonPoints
{
    class First50PrimesGenerator : IPrimeGenerator
    {
        private static readonly int[] primes = 
          { 2,   3,   5,   7,   11,  13,  17,  19,  23,  29,
            31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
            73,  79,  83,  89,  97,  101, 103, 107, 109, 113,
            127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
            179, 181, 191, 193, 197, 199, 211, 223, 227, 229};

        public int[] FirstPrimes(int count)
        {
            if (count > 50)
            {
                throw new ArgumentException("This sieve can generate only the first 50 primes");
            }
            return (int[])primes.Clone();
        }
    }
}