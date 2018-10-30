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

        public RandomVariable(double magnitude)
        {
            Magnitude = magnitude;
        }

        public double Realize(int iteration, int parameters)
        {
            return random.NextDouble() * Magnitude ;
        }
    }
}
