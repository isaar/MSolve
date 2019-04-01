using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using Troschuetz.Random;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Analyzers
{
    public enum RandomVariableDistributionType
    {
        Normal,
        Lognormal
    }

    public class RandomVariableTargetEvaluator:IStochasticMaterialCoefficientsProvider
    {
        private readonly double _standardDeviation;
        private readonly double _meanValue;
        private double _randomVariable;
        private readonly RandomVariableDistributionType _distributionType;


        /// <summary>Initializes a new instance of the <see cref="RandomVariableTargetEvaluator"/> class.</summary>
        /// <param name="meanValue">The mean value.</param>
        /// <param name="standardDeviation">The standard deviation.</param>
        /// <param name="distributionType">Type of the distribution.</param>
        public RandomVariableTargetEvaluator(double meanValue, double standardDeviation, RandomVariableDistributionType distributionType)
        {
            this._meanValue = meanValue;
            this._standardDeviation = standardDeviation;
            _distributionType = distributionType;
        }

        /// <summary>Gets or sets the random variables.</summary>
        /// <value>The random variables.</value>
        public double[] RandomVariables
        {
            get { return new double[1];}
            set { CalculateRandomVariable(); }
        }

        /// <summary>
        ///   <para>
        ///  Calculates the random variable.
        /// </para>
        /// </summary>
        private void CalculateRandomVariable()
        {
            switch (_distributionType)
            {
                case RandomVariableDistributionType.Normal:
                    var normalDistribution = new NormalDistribution(_meanValue, _standardDeviation);
                    _randomVariable = normalDistribution.NextDouble();
                    break;
                case RandomVariableDistributionType.Lognormal:
                    var lognormal = new LognormalDistribution(_meanValue, _standardDeviation*_meanValue);
                    _randomVariable = lognormal.NextDouble();
                    break;
            }
        }

        /// <summary>Gets the coefficient.</summary>
        /// <param name="meanValue">The mean value.</param>
        /// <param name="coordinates">The coordinates.</param>
        /// <returns></returns>
        public double GetCoefficient(double meanValue, double[] coordinates)
        {
            return 1/_randomVariable;
        }
    }
}
