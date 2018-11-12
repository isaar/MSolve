using MGroup.Stochastic.Interfaces;
using System;

namespace MGroup.Stochastic
{
    public class MonteCarlo
    {
        public int NoOfIterations { get; }
        public ISystemRealizer SystemRealizer { get; }
        public ISystemResponseEvaluator SystemResponseEvaluator { get; }

        public MonteCarlo(int noOfIterations, ISystemRealizer systemRealizer, ISystemResponseEvaluator systemResponseEvaluator)
        {
            NoOfIterations = noOfIterations;
            SystemRealizer = systemRealizer;
            SystemResponseEvaluator = systemResponseEvaluator;
        }

        public void Evaluate()
        {
            SystemRealizer.Realize(0);
            int systemResponseDimension = SystemResponseEvaluator.Evaluate(0).Length;
            var systemResponse = new double[systemResponseDimension];
            var MonteCarloMeanValue = new double[systemResponseDimension];
            var MonteCarloStandardDeviation = new double[systemResponseDimension];
            systemResponse = SystemResponseEvaluator.Evaluate(0);
            for (int i = 0; i < systemResponse.Length; i++)
            {
                MonteCarloMeanValue[i] += systemResponse[i];
            }

            for (int i = 1; i < NoOfIterations; i++)
            {
                SystemRealizer.Realize(i);
                systemResponse = SystemResponseEvaluator.Evaluate(i);
                for (int j = 0; j < systemResponse.Length; j++)
                {
                    MonteCarloMeanValue[j] += systemResponse[j];
                }
            }

            for (int i = 0; i < systemResponse.Length; i++)
            {
                MonteCarloMeanValue[i] = MonteCarloMeanValue[i] / NoOfIterations;
                MonteCarloStandardDeviation[i] = (systemResponse[i] - MonteCarloMeanValue[i]) * (systemResponse[i] - MonteCarloMeanValue[i]);
                MonteCarloStandardDeviation[i] = Math.Sqrt(MonteCarloStandardDeviation[i] / NoOfIterations);
            }
        }
    }
}
