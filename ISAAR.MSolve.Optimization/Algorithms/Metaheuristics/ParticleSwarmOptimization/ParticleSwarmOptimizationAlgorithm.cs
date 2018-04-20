using ISAAR.MSolve.Optimization.Convergence;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Problems;
using System;
using static ISAAR.MSolve.Optimization.Commons.VectorOperations;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.ParticleSwarmOptimization
{
    /// <summary>
    ///     Class implementing the Particle Swarm Optimization metaheuristic search algorithm
    ///     <see href = " https://doi.org/10.1109/ICNN.1995.488968"> 
    ///      J. Kennedy and R. Eberhart (1995). "Particle Swarm Optimization". 
    ///     Proceedings of IEEE International Conference on Neural Networks. pp. 1942–1948.
    ///     </see>
    /// </summary>
    public class ParticleSwarmOptimizationAlgorithm : IOptimizationAlgorithm
    {
        private readonly int swarmSize;
        private readonly int dimension;

        private readonly double omega;                 // Inertia weight
        private readonly double phip;                  // Scaling factor to search away from the particle's best known position
        private readonly double phig;                  // Scaling factor to search away from the swarm's best known position

        private readonly IDesignFactory designFactory;
        private readonly double[] lowerBound;
        private readonly double[] upperBound;
        private readonly IConvergenceCriterion convergenceCriterion;
        private readonly IOptimizationLogger logger;

        private Particle[] individuals;

        private readonly Random randomNumberGenerator = new Random();

        public ParticleSwarmOptimizationAlgorithm(int dimension, double[] lowerBound, double[] upperBound, IDesignFactory designFactory,
           int swarmSize, double omega, double phip, double phig, IConvergenceCriterion convergenceCriterion, IOptimizationLogger logger)
        {
            this.dimension = dimension;
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
            this.designFactory = designFactory;

            this.swarmSize = swarmSize;
            this.omega = omega;
            this.phip = phip;
            this.phig = phig;
            this.convergenceCriterion = convergenceCriterion;
            this.logger = logger;
        }

        public double BestFitness
        {
            get; private set;
        } = double.MaxValue;

        public double[] BestPosition
        {
            get; private set;
        }


        public int CurrentIteration
        {
            get; private set;
        }

        public double CurrentFunctionEvaluations
        {
            get; private set;
        }

        private void Initialize()
        {
            individuals = new Particle[swarmSize];
            double[] vhigh = new double[dimension];
            double[] vlow = new double[dimension];

            // Initialize the swarm's position with a uniformly distributed random vector as well as the velocity.
            for (int i = 0; i < swarmSize; i++)
            {
                double[] position = new double[dimension];
                double[] velocity = new double[dimension];
                double[] bestPosition = new double[dimension];
                double[] personalBestPosition = new double[dimension];

                for (int v = 0; v < dimension; v++)
                {
                    vhigh[v] = Math.Abs(upperBound[v] - lowerBound[v]);
                    vlow[v] = -vhigh[v];

                    position[v] = lowerBound[v] + randomNumberGenerator.NextDouble() * (upperBound[v] - lowerBound[v]);
                    velocity[v] = vlow[v] + randomNumberGenerator.NextDouble() * (vhigh[v] - vlow[v]);
                    bestPosition[v] = position[v];
                }
                individuals[i] = new Particle(position, velocity, double.MaxValue, personalBestPosition, double.MaxValue);
            }

            // Evaluate the initial population
            Evaluation(individuals);

            // Store the best fitness and position
            for (int i = 0; i < swarmSize; i++)
            {
                if (individuals[i].ObjectiveValue < BestFitness)
                {
                    BestFitness = individuals[i].ObjectiveValue;
                    BestPosition = individuals[i].Position;
                }
            }
        }

        private void Evaluation(Individual[] individuals)
        {
            for (int i = 0; i < swarmSize; i++)
            {
                IDesign design = designFactory.CreateDesign(individuals[i].Position);
                double fitness = design.ObjectiveValues[0];

                individuals[i].ObjectiveValue = fitness;
            }
            CurrentFunctionEvaluations += swarmSize;
        }

        private void Iterate()
        {
            UpdatePositionAndVelocity();
            CheckBounds();
            Evaluation(individuals);
            UpdateBest();
        }

        private void UpdatePositionAndVelocity()
        {

            for (int i = 0; i < swarmSize; i++)
            {
                int rp = randomNumberGenerator.Next(dimension);
                int rg = randomNumberGenerator.Next(dimension);

                individuals[i].Velocity = Add(Add(Scale(omega, individuals[i].Velocity), 
                                                  Scale(phip, Scale(rp, Subtract(BestPosition, individuals[i].Position)))), 
                                              Scale(phig, Scale(rg, Subtract(individuals[i].PersonalBestPosition, individuals[i].Position))));

                individuals[i].Position = Add(individuals[i].Position, individuals[i].Velocity);
            }
        }

        /// <summary>
        ///     Verify upper and lower bounds for the updated position vector
        /// </summary>
        private void CheckBounds()
        {
            for (int i = 0; i < swarmSize; i++)
            {
                double[] PositionVector = individuals[i].Position;

                for (int j = 0; j < dimension; j++)
                {
                    if (PositionVector[j] > upperBound[j]) PositionVector[j] = upperBound[j];
                    else if (PositionVector[j] < lowerBound[j]) PositionVector[j] = lowerBound[j];
                }
            }
        }

        private void UpdateBest()
        {
            for (int i = 0; i < swarmSize; i++)
            {
                if (individuals[i].ObjectiveValue < individuals[i].PersonalBestFitness)
                {
                    individuals[i].PersonalBestPosition = individuals[i].Position;
                    individuals[i].PersonalBestFitness = individuals[i].ObjectiveValue;
                }

                if (individuals[i].ObjectiveValue < BestFitness)
                {
                    BestPosition = individuals[i].Position;
                    BestFitness = individuals[i].ObjectiveValue;
                }
            }
        }

        public void Solve()
        {
            Initialize();
            logger.Log(this);

            CurrentIteration = 0;
            while (!convergenceCriterion.HasConverged(this))
            {
                CurrentIteration++;
                Iterate();
                logger.Log(this);

                // Write best fitness and position
                Console.WriteLine(String.Format(@"Iter: {0} | {1} ", CurrentIteration, BestFitness));
                // Array.ForEach(BestPosition, x => Console.WriteLine(x));
            }
        }

        public class Builder
        {
            public int SwarmSize { get; set; }
            public double Omega { get; set; } = 0.5;
            public double PhiP { get; set; } = 2.0;
            public double PhiG { get; set; } = 2.0;
            public IConvergenceCriterion ConvergenceCriterion { get; set; }
            public IOptimizationLogger Logger { get; set; }

            private OptimizationProblem optimizationProblem;

            public Builder(OptimizationProblem optimizationProblem)
            {
                ProblemChecker.Check(optimizationProblem);
                this.optimizationProblem = optimizationProblem;
            }

            public IOptimizationAlgorithm Build()
            {
                return new ParticleSwarmOptimizationAlgorithm(optimizationProblem.Dimension,
                    optimizationProblem.LowerBound, optimizationProblem.UpperBound, optimizationProblem.DesignFactory,
                    SwarmSize, Omega, PhiP, PhiG, ConvergenceCriterion, Logger);
            }
        }
    }
}