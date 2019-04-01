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

        /// <summary>Initializes a new instance of the <see cref="RandomVariable"/> class with specified mean value, standard deviation and pdf.</summary>
        /// <param name="meanValue">The mean value.</param>
        /// <param name="standardDeviation">The standard deviation.</param>
        /// <param name="distributionType">Type of the distribution.</param>
        public RandomVariable(double meanValue, double standardDeviation, RandomVariableDistributionType distributionType)
        {
            this._meanValue = meanValue;
            this._standardDeviation = standardDeviation;
            _distributionType = distributionType;
        }

        /// <summary>Realizes the specified iteration for given stochasti domain mapper and domain prameters.</summary>
        /// <param name="iteration">The iteration.</param>
        /// <param name="domainMapper">The domain mapper.</param>
        /// <param name="parameters">The parameters.</param>
        /// <returns></returns>
        /// <exception cref="System.NotImplementedException"></exception>
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
