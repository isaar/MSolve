using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Initialization.HaltonPoints;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Initialization
{
    public class RealUniformHaltonInitializer: IInitializer<double>
    {
        private readonly HaltonPointGenerator pointGenerator;
        private readonly int continuousVariablesCount;
        private readonly double[] lowerBounds;
        private readonly double[] upperBounds;

        public RealUniformHaltonInitializer(OptimizationProblem problem): 
                            this(problem, new HaltonPointGenerator(problem.Dimension))
        {
        }

        public RealUniformHaltonInitializer(OptimizationProblem problem, 
            HaltonPointGenerator pointGenerator)
        {
            this.continuousVariablesCount = problem.Dimension;
            this.lowerBounds = problem.LowerBound;
            this.upperBounds = problem.UpperBound;

            if (pointGenerator == null)
            {
                throw new ArgumentException("The Halton point generator must not be null");
            }
            this.pointGenerator = pointGenerator;
        }

        public double[] Generate()
        {
            double[] sample = pointGenerator.NextPoint();
            for (int i = 0; i < continuousVariablesCount; ++i)
            {
                // min+rand*(max-min) could overflow. This won't
                sample[i] = (1 - sample[i]) * lowerBounds[i] + sample[i] * upperBounds[i]; 
            }
            return sample;
        }
    }
}
