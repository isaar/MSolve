using MGroup.Stochastic.Interfaces;
using System;

namespace MGroup.Stochastic
{
    public class MonteCarlo
    {
        public int NoOfIterations { get; }
        public double[] MeanValueResponse;
        public double[] StandardDeviationResponse;
        public ISystemRealizer SystemRealizer { get; }
        public ISystemResponseEvaluator SystemResponseEvaluator { get; }

        /// <summary>Initializes a new instance of the <see cref="MonteCarlo"/> class.</summary>
        /// <param name="noOfIterations">The no of iterations.</param>
        /// <param name="systemRealizer">The system realizer.</param>
        /// <param name="systemResponseEvaluator">The system response evaluator.</param>
        public MonteCarlo(int noOfIterations, ISystemRealizer systemRealizer, ISystemResponseEvaluator systemResponseEvaluator)
        {
            NoOfIterations = noOfIterations;
            SystemRealizer = systemRealizer;
            SystemResponseEvaluator = systemResponseEvaluator;
        }

        /// <summary>Evaluates 1st and 2nd statistical moments of the predesignated response.</summary>
        public void Evaluate()
        {
            SystemRealizer.Realize(0);
            int systemResponseDimension = SystemResponseEvaluator.Evaluate(0).Length;
            var systemResponse = new double[systemResponseDimension];
            MeanValueResponse = new double[systemResponseDimension];
            StandardDeviationResponse = new double[systemResponseDimension];
            systemResponse = SystemResponseEvaluator.Evaluate(0);
            for (int i = 0; i < systemResponse.Length; i++)
            {
                MeanValueResponse[i] += systemResponse[i];
            }

            for (int i = 1; i < NoOfIterations; i++)
            {
                SystemRealizer.Realize(i);
                systemResponse = SystemResponseEvaluator.Evaluate(i);
                for (int j = 0; j < systemResponse.Length; j++)
                {
                    MeanValueResponse[j] += systemResponse[j];
                }
            }

            for (int i = 0; i < systemResponse.Length; i++)
            {
                MeanValueResponse[i] = MeanValueResponse[i] / NoOfIterations;
                StandardDeviationResponse[i] = (systemResponse[i] - MeanValueResponse[i]) * (systemResponse[i] - MeanValueResponse[i]);
                StandardDeviationResponse[i] = Math.Sqrt(StandardDeviationResponse[i] / NoOfIterations);
            }
        }
    }
}
