using MGroup.Stochastic.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.Stochastic.Structural.StochasticRealizers
{
    public class RandomVariable : IUncertainParameterRealizer
    {
        private Random random = new Random();
        public double Magnitude { get; }
        public IStochasticDomainMapper DomainMapper;
        public RandomVariable(double magnitude, IStochasticDomainMapper domainMapper)
        {
            Magnitude = magnitude;
            DomainMapper = domainMapper;
        }

        public double Realize(int iteration, IStochasticDomainMapper domainMapper, double[] parameters)
        {
            return random.NextDouble() * Magnitude ;
        }
    }
}
