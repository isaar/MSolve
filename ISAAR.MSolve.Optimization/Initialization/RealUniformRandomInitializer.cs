using ISAAR.MSolve.Optimization.Commons;
using ISAAR.MSolve.Optimization.Problems;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Initialization
{
    public class RealUniformRandomInitializer : IInitializer<double>
    {
        private readonly IGenerator rng;
        private readonly int continuousVariablesCount;
        private readonly double[] lowerBounds;
        private readonly double[] upperBounds;

        public RealUniformRandomInitializer(OptimizationProblem problem) : 
                            this(problem, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public RealUniformRandomInitializer(OptimizationProblem problem, IGenerator randomNumberGenerator)
        {
            //problem.CheckInput();
            this.continuousVariablesCount = problem.Dimension;
            this.lowerBounds = problem.LowerBound;
            this.upperBounds = problem.UpperBound;

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.rng = randomNumberGenerator;
        }

        public double[] Generate()
        {
            double[] sample = new double[continuousVariablesCount];
            for (int i = 0; i < continuousVariablesCount; ++i)
            {
                double rand = rng.NextDouble();
                sample[i] = (1 - rand) * lowerBounds[i] + rand * upperBounds[i]; // min+rand*(max-min) could overflow. This won't
            }
            return sample;
        }
    }
}
