using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings
{
    internal class Quantization
    {
        private readonly int bitsCount;
        private readonly double[] subranges;

        internal Quantization(int bitsCount)
        {
            this.bitsCount = bitsCount;
            int decimalsCount = (int)Math.Pow(2, bitsCount);
            subranges = new double[decimalsCount + 1];
            for (int i = 0; i < subranges.Length; ++i)
            {
                subranges[i] = i * (1.0 / (decimalsCount));
            }
            subranges[decimalsCount] = 1.0;
        }

        internal int NormalizedDoubleToInteger(double normalizedDouble)
        {
            if (normalizedDouble < 0 || normalizedDouble >= 1)
            {
                throw new ArgumentException("The provided value: " + normalizedDouble + " is not normalized");
            }
            for (int i = 1; i < subranges.Length; ++i)
            {
                if (normalizedDouble < subranges[i]) return i - 1;
            }
            throw new Exception("This should not have been reached");
        }

        internal void PrintRanges()
        {
            Console.WriteLine("Subranges: ");
            for (int i = 0; i < subranges.Length-1; ++i)
            {
                Console.WriteLine(subranges[i] + " - " + subranges[i+1]);
            }
        }


    }
}
