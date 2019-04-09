using System;
using MGroup.Stochastic.Interfaces;

namespace ISAAR.MSolve.Tests.SupportiveClasses
{
    public class CantileverStochasticDomainMapper : IStochasticDomainMapper
    {
        private readonly double[] origin;

        public CantileverStochasticDomainMapper(double[] origin)
        {
            this.origin = origin;
        }

        public double[] Map(double[] problemDomainVector)
        {
            return new[]
            {
                Math.Sqrt(Math.Pow(problemDomainVector[0] - origin[0], 2) +
                          Math.Pow(problemDomainVector[1] - origin[1], 2) +
                          Math.Pow(problemDomainVector[2] - origin[2], 2))
            };
        }
    }
}
