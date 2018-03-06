using ISAAR.MSolve.Optimization.Commons;
using ISAAR.MSolve.Optimization.Problems;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations.Gaussian
{
    public class ConstantGaussianMutation : AbstractPerturbationTemplate
    {
        private readonly int continuousVariablesCount;
        private readonly double[] standardDeviations;

        // If the user needs reproductible gaussians he will need to pass in a "rigged" IGenerator 
        // and the NormalDistribution will be composed with that IGenerator
        private readonly NormalDistribution normalDistribution;

        public ConstantGaussianMutation(OptimizationProblem problem, double standardDeviation = 1.0):
                    this(problem, standardDeviation, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public ConstantGaussianMutation(OptimizationProblem problem, double[] standardDeviations) :
                    this(problem, standardDeviations, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public ConstantGaussianMutation(OptimizationProblem problem, double standardDeviation, IGenerator randomNumberGenerator)
        {
            //problem.CheckInput();
            this.continuousVariablesCount = problem.Dimension;

            if (standardDeviation <= 0)
            {
                throw new ArgumentException("Standard deviation must be > 0, but was " + standardDeviation);
            }
            this.standardDeviations = new double[continuousVariablesCount];
            for (int i = 0; i < problem.Dimension; ++i)
            {
                this.standardDeviations[i] = standardDeviation;
            }

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.normalDistribution = new NormalDistribution(randomNumberGenerator, 0, 1);
        }

        public ConstantGaussianMutation(OptimizationProblem problem, double[] standardDeviations, IGenerator randomNumberGenerator)
        {
            //problem.CheckInput();
            this.continuousVariablesCount = problem.Dimension;

            if (standardDeviations.Length != continuousVariablesCount)
            {
                throw new ArgumentException("There are " + continuousVariablesCount + " continuous design variables, but "
                                            + standardDeviations.Length + " standard deviations are provided");
            }
            this.standardDeviations = new double[continuousVariablesCount];
            for (int i = 0; i < problem.Dimension; ++i)
            {
                if (standardDeviations[i] <= 0)
                {
                    throw new ArgumentException("Standard deviations must be > 0, but at index " + i 
                                                + ", the provided value is " + standardDeviations[i]);
                }
                this.standardDeviations[i] = standardDeviations[i];
            }

            if (randomNumberGenerator == null) throw new ArgumentException("The random number generator must not be null");
            this.normalDistribution = new NormalDistribution(randomNumberGenerator, 0, 1);
        }

        protected override double[] ComputePerturbations()
        {
            double[] perturbations = new double[continuousVariablesCount];
            for (int i = 0; i < continuousVariablesCount; ++i)
            {
                perturbations[i] = standardDeviations[i] * normalDistribution.NextDouble();
            }
            return perturbations;
        }
    }
}
