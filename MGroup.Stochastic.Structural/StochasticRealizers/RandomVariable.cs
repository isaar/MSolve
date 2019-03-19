using MGroup.Stochastic.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using PommaLabs.Thrower;
using Troschuetz.Random;
using Troschuetz.Random.Distributions.Continuous;

namespace MGroup.Stochastic.Structural.StochasticRealizers
{
    public enum RandomVariableDistributionType
    {
        Normal,
        Lognormal
    }
    public class RandomVariable : IUncertainParameterRealizer
    {
        private readonly double _standardDeviation;
        private readonly double _meanValue;
        private double _randomVariable;
        private readonly RandomVariableDistributionType _distributionType;
        public IStochasticDomainMapper DomainMapper;

        public RandomVariable(double meanValue, double standardDeviation, RandomVariableDistributionType distributionType)
        {
            this._meanValue = meanValue;
            this._standardDeviation = standardDeviation;
            _distributionType = distributionType;
        }

        public double Realize(int iteration, IStochasticDomainMapper domainMapper, double[] parameters)
        {
            switch (_distributionType)
            {
                case RandomVariableDistributionType.Normal:
                    var normalDistribution = new NormalDistribution(_meanValue, _standardDeviation);
                    return _randomVariable = normalDistribution.NextDouble();
                    break;
                case RandomVariableDistributionType.Lognormal:
                    var lognormal = new LognormalDistribution(_meanValue, _standardDeviation * _meanValue);
                    return _randomVariable = lognormal.NextDouble();
                    break;
                default:
                    return 0;
                    throw new NotImplementedException();
                    break;
            }
        }

    }
}
