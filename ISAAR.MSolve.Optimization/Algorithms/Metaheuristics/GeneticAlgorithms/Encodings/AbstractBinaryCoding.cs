using ISAAR.MSolve.Optimization.Commons;
using ISAAR.MSolve.Optimization.Problems;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings
{
    public abstract class AbstractBinaryCoding: IEncoding<bool>
    {
        #region fields and properties
        private readonly IGenerator rng = RandomNumberGenerationUtilities.troschuetzRandom;
        private readonly Quantization quantization;

        // Continuous variables
        private readonly int continuousVariablesCount;
        private readonly double[] continuousLowerBounds;
        private readonly double[] continuousUpperBounds;
        private readonly int bitsPerContinuousVariable;

        // Integer variables
        private readonly int integerVariablesCount;
        //private readonly int[] integerLowerBounds; They are 0
        //private readonly int[] integerUpperBounds;
        private readonly int bitsPerIntegerVariable;
        #endregion

        #region constructor
        protected AbstractBinaryCoding(OptimizationProblem problem, int bitsPerContinuousVariable, int bitsPerIntegerVariable)
        {
            this.continuousVariablesCount = problem.Dimension;
            this.continuousLowerBounds = problem.LowerBound;
            this.continuousUpperBounds = problem.UpperBound;
            this.bitsPerContinuousVariable = bitsPerContinuousVariable;
            this.integerVariablesCount = 0;
            //this.integerUpperBounds = null;
            this.bitsPerIntegerVariable = 0;

            this.quantization = new Quantization(bitsPerContinuousVariable);
        }
        #endregion

        #region IEncoding implementations
        public bool[] ComputeGenotype(double[] phenotype)
        {
            // Continuous variables
            bool[] genotype = new bool[continuousVariablesCount * bitsPerContinuousVariable];
            for (int i = 0; i < continuousVariablesCount; ++i)
            {
                double normalized = (phenotype[i] - continuousLowerBounds[i]) / (continuousUpperBounds[i] - continuousLowerBounds[i]);
                int dec = quantization.NormalizedDoubleToInteger(normalized);
                DecimalIntegerToBitstring(dec, genotype, i * bitsPerContinuousVariable, bitsPerContinuousVariable);
            }
            return genotype;
        }

        public double[] ComputePhenotype(bool[] genotype)
        {
            // Continuous variables
            double[] continuousVariables = new double[continuousVariablesCount];
            for (int i = 0; i < continuousVariablesCount; ++i)
            {
                int start = i * bitsPerContinuousVariable;
                int dec = BitstringToDecimalInteger(genotype, start, bitsPerContinuousVariable);
                double normalized = dec / (Math.Round(Math.Pow(2, bitsPerContinuousVariable)) - 1);
                continuousVariables[i] = continuousLowerBounds[i] +
                                         normalized * (continuousUpperBounds[i] - continuousLowerBounds[i]);
            }
            return continuousVariables;
        }

        public int[] IntegerPhenotype(bool[] genotype)
        {
            // Integer Variables
            int[] integerVariables = new int[integerVariablesCount];
            int offset = continuousVariablesCount * bitsPerContinuousVariable;
            for (int i = 0; i < integerVariablesCount; ++i)
            {
                int start = offset + i * bitsPerIntegerVariable;
                integerVariables[i] = BitstringToDecimalInteger(genotype, start, bitsPerIntegerVariable);
            }
            return integerVariables;
        }
        #endregion

        #region abstract methods
        protected abstract int BitstringToDecimalInteger(bool[] bits, int start, int length);
        protected abstract void DecimalIntegerToBitstring(int dec, bool[] bits, int start, int length);

        #endregion
    }
}
