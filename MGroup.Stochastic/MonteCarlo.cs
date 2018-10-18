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
            for (int i = 0; i < NoOfIterations; i++)
            {
                SystemRealizer.Realize(i);
                SystemResponseEvaluator.Evaluate(i);
            }
        }
    }
}
