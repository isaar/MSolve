namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.ParticleSwarmOptimization
{
    public class Particle : Individual
    {

        public double[] Velocity
        {
            get; set;
        }

        public double[] PersonalBestPosition
        {
            get; set;
        }

        public double PersonalBestFitness
        {
            get; set;
        } = double.MaxValue;

        public Particle(double[] position, double[] velocity, double objectiveValue, 
            double[] personalBestPosition, double personalBestFitness) : base(position, objectiveValue)
        {
            this.PersonalBestPosition = personalBestPosition;
            this.PersonalBestFitness = personalBestFitness;
            this.Velocity = velocity;
        }
    }
}